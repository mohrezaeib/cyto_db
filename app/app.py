from flask import Flask, render_template, request, abort, url_for
import math
import json

app = Flask(__name__)


# Load data from JSON file
def load_data():
    with open("data.json", "r", encoding="utf-8") as file:
        return json.load(file)

data = load_data()

ITEMS_PER_PAGE = 10

ITEMS_PER_PAGE = 10

@app.route("/")
def list_items():
    # Get search query and selected fields from the URL parameters.
    query = request.args.get("query", "").strip()
    selected_fields = request.args.getlist("fields")  # may be empty

    # Get current page (default 1)
    try:
        page = int(request.args.get("page", 1))
    except ValueError:
        page = 1

    # Filtering: if a search query exists, filter based on the selected fields.
    filtered_data = data
    if query:
        filtered_data = []
        for item in data:
            # If specific fields are selected, only search in those;
            # otherwise search in all fields.
            fields_to_search = selected_fields or item["fields"].keys()
            for field in fields_to_search:
                value = item["fields"].get(field, "")
                if query.lower() in str(value).lower():
                    filtered_data.append(item)
                    break  # if one field matches, no need to check others

    # Pagination logic.
    total_items = len(filtered_data)
    total_pages = math.ceil(total_items / ITEMS_PER_PAGE)
    start = (page - 1) * ITEMS_PER_PAGE
    end = start + ITEMS_PER_PAGE
    paginated_data = filtered_data[start:end]

    # Use keys from the first item's fields as possible search options.
    possible_fields = list(data[0]["fields"].keys()) if data else []

    return render_template("list.html",
                           items=paginated_data,
                           page=page,
                           total_pages=total_pages,
                           query=query,
                           selected_fields=selected_fields,
                           possible_fields=possible_fields)

@app.route("/item/<int:mol_idx>")
def item_detail(mol_idx):
    # Find the item by mol_idx.
    item = next((item for item in data if item["mol_idx"] == mol_idx), None)
    if not item:
        abort(404)
    return render_template("detail.html", item=item)

if __name__ == "__main__":
    # Expose the app on all network interfaces at port 5000.
    app.run(host="0.0.0.0", port=5000, debug=True)
