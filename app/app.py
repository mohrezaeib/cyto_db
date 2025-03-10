from flask import Flask, render_template, request, jsonify
import json
from config  import ALLOWED_FIELDS , ALLOWED_LIST_FIELDS # Import allowed fields

app = Flask(__name__)

# Load data from a local JSON file
with open("data.json", "r") as f:
    data = json.load(f)

# Pagination settings
ITEMS_PER_PAGE = 10

@app.route("/")
def list_items():
    query = request.args.get("query", "").strip()
    selected_fields = request.args.getlist("fields")

    # Default to all allowed fields if none are selected
    if not selected_fields:
        selected_fields = ALLOWED_FIELDS

    # Filter items based on search
    filtered_data = []
    for item in data:
        if not query:
            filtered_data.append(item)
        else:
            for field in selected_fields:
                if field in item["fields"] and query.lower() in str(item["fields"][field]).lower():
                    filtered_data.append(item)
                    break  # Avoid duplicates

    # Pagination logic
    page = int(request.args.get("page", 1))
    total_pages = (len(filtered_data) + ITEMS_PER_PAGE - 1) // ITEMS_PER_PAGE
    paginated_items = filtered_data[(page - 1) * ITEMS_PER_PAGE: page * ITEMS_PER_PAGE]

    return render_template(
        "list.html",
        items=paginated_items,
        page=page,
        total_pages=total_pages,
        query=query,
        all_list_fields=ALLOWED_LIST_FIELDS,  # Pass allowed fields to template
        all_fields=ALLOWED_FIELDS,  # Pass allowed fields to template
        selected_fields=selected_fields,
    )

@app.route("/item/<int:mol_idx>")
def item_detail(mol_idx):
    # Find current item
    item = next((i for i in data if i["mol_idx"] == mol_idx), None)
    if not item:
        return "Item not found", 404

    # Find index of the current item
    idx = next((i for i, d in enumerate(data) if d["mol_idx"] == mol_idx), None)

    # Determine previous and next items
    prev_idx = data[idx - 1]["mol_idx"] if idx   > 0 else None
    next_idx = data[idx + 1]["mol_idx"] if idx < len(data) - 1 else None

    return render_template("detail.html", item=item, prev_idx=prev_idx, next_idx=next_idx)

