from flask import Flask, render_template, request, abort, url_for
import json
import math

app = Flask(__name__)

# Load data from JSON file
def load_data():
    with open("data.json", "r", encoding="utf-8") as file:
        return json.load(file)

data = load_data()

ITEMS_PER_PAGE = 10

def get_all_fields(data):
    """Return a sorted list of all unique field keys from the data."""
    all_fields = set()
    for item in data:
        all_fields.update(item.get("fields", {}).keys())
    return sorted(all_fields)

@app.route("/")
def list_items():
    query = request.args.get("query", "").strip()
    selected_fields = request.args.getlist("fields")

    try:
        page = int(request.args.get("page", 1))
    except ValueError:
        page = 1

    filtered_data = data
    if query:
        filtered_data = []
        for item in data:
            # If no fields are selected, search in all available fields
            fields_to_search = selected_fields or item["fields"].keys()
            for field in fields_to_search:
                value = item["fields"].get(field, "")
                if query.lower() in str(value).lower():
                    filtered_data.append(item)
                    break

    total_items = len(filtered_data)
    total_pages = math.ceil(total_items / ITEMS_PER_PAGE)
    start = (page - 1) * ITEMS_PER_PAGE
    end = start + ITEMS_PER_PAGE
    paginated_data = filtered_data[start:end]

    all_fields = get_all_fields(data)

    return render_template("list.html",
                           items=paginated_data,
                           page=page,
                           total_pages=total_pages,
                           query=query,
                           selected_fields=selected_fields,
                           all_fields=all_fields)

@app.route("/item/<int:mol_idx>")
def item_detail(mol_idx):
    item = next((item for item in data if item["mol_idx"] == mol_idx), None)
    if not item:
        abort(404)
    return render_template("detail.html", item=item)

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000, debug=True)
