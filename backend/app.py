from flask import Flask, request, jsonify
from flask_cors import CORS
from flask_caching import Cache
import json
import os
import logging
import pandas as pd
import numpy as np
from dotenv import load_dotenv

# Ensure these match your project structure
# from params import FilterParams  <-- Assuming this exists, but defining a mock below for safety if you copy-paste
from dataclasses import dataclass, field
from typing import List, Optional

load_dotenv()

# --- CONFIGURATION ---
BACKEND_PORT = int(os.getenv('CYTOCHAL_BACKEND_PORT', 5000))
CYTOCHAL_API_BASE_PATH = os.getenv("CYTOCHAL_API_BASE_PATH", "cytochal-api")
FLASK_ENV = os.getenv('CYTOCHAL_FLASK_ENV', 'production')
DEBUG = os.getenv('CYTOCHAL_DEBUG', 'False').lower() == 'true'
DATA_PATH = os.getenv('CYTOCHAL_DATA_PATH', 'static/data.json')
LOG_LEVEL = os.getenv('CYTOCHAL_LOG_LEVEL', 'INFO')

# --- LOGGING ---
logging.basicConfig(level=getattr(logging, LOG_LEVEL.upper(), logging.INFO))
logger = logging.getLogger(__name__)

app = Flask(__name__)
CORS(app)

# --- CACHE ---
app.config['CACHE_TYPE'] = 'SimpleCache'
app.config['CACHE_DEFAULT_TIMEOUT'] = 300
cache = Cache(app)

# --- DATACLASS (Mocking your existing params.py for completeness) ---
@dataclass
class FilterParams:
    query: str = ""
    page: int = 1
    per_page: int = 20
    min_molweight: Optional[float] = None
    max_molweight: Optional[float] = None
    min_ic50: Optional[float] = None
    max_ic50: Optional[float] = None
    activity: Optional[str] = None
    reversibility: Optional[str] = None
    quantity_type: Optional[str] = None
    min_quantity: Optional[float] = None
    max_quantity: Optional[float] = None
    selected_fields: List[str] = field(default_factory=list)

# --- PANDAS DATA LOADING ---
# Column Mapping: Internal Name -> JSON Field Key (flattened)
COL_MAP = {
    "weight": "fields.Total Molweight",
    "ic50": "fields.Cytotoxicity in L929 mouse fibroblasts [µM]",
    "activity": "fields.Actin Disruption Activity (1h Treatment)",
    "reversibility": "fields.Reversibilty", # Note: Preserving source typo
    "quantity": "fields.Quantity",
    "compound": "fields.Compound",
    "smiles": "fields.Smiles"
}

def load_and_prep_data():
    if not os.path.exists(DATA_PATH):
        raise FileNotFoundError(f"Data file not found at: {DATA_PATH}")

    with open(DATA_PATH, "r", encoding="utf-8") as f:
        raw_data = json.load(f)

    # 1. Flatten JSON: propels 'fields' up to top level with dot notation
    df = pd.json_normalize(raw_data)
    
    # 2. Sort by mol_idx to ensure consistent Next/Prev navigation
    if 'mol_idx' in df.columns:
        df = df.sort_values('mol_idx').reset_index(drop=True)

    # 3. Clean Numeric Columns (Coerce errors turns "not tested" into NaN)
    #    This prevents crashes when filtering strings as numbers
    for col in [COL_MAP['weight'], COL_MAP['ic50']]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    
    logger.info(f"Successfully loaded and processed {len(df)} compounds.")
    return df

# Initialize Global DataFrame
try:
    DF = load_and_prep_data()
except Exception as e:
    logger.error(f"Failed to load data: {e}")
    DF = pd.DataFrame() # Fallback empty DF

# --- ROUTES ---

@app.route("/cytochal-api/health")
def health_check():
    return jsonify({
        "status": "healthy", 
        "data_loaded": len(DF),
        "environment": FLASK_ENV
    })

@app.route(f"/{CYTOCHAL_API_BASE_PATH}/items")
def api_items():
    try:
        # Parse Params
        params = FilterParams(
            query=request.args.get("query", "").strip().lower(),
            page=max(1, request.args.get("page", default=1, type=int)),
            per_page=request.args.get("per_page", default=20, type=int),
            min_molweight=request.args.get("min_molweight", type=float),
            max_molweight=request.args.get("max_molweight", type=float),
            min_ic50=request.args.get("min_ic50", type=float),
            max_ic50=request.args.get("max_ic50", type=float),
            activity=request.args.get("activity"),
            reversibility=request.args.get("reversibility"),
            quantity_type=request.args.get("quantity_type"),
            min_quantity=request.args.get("min_quantity", type=float),
            max_quantity=request.args.get("max_quantity", type=float),
            selected_fields=request.args.getlist("fields")
        )

        # Start with a mask of all True
        mask = pd.Series(True, index=DF.index)

        # --- A. Numeric Range Filters (Vectorized) ---
        if params.min_molweight: mask &= (DF[COL_MAP['weight']] >= params.min_molweight)
        if params.max_molweight: mask &= (DF[COL_MAP['weight']] <= params.max_molweight)
        if params.min_ic50:      mask &= (DF[COL_MAP['ic50']] >= params.min_ic50)
        if params.max_ic50:      mask &= (DF[COL_MAP['ic50']] <= params.max_ic50)

        # --- B. Categorical Filters (Exact Match) ---
        if params.activity:      mask &= (DF[COL_MAP['activity']] == params.activity)
        if params.reversibility: mask &= (DF[COL_MAP['reversibility']] == params.reversibility)

        # --- C. Text Search ---
        if params.query:
            if params.selected_fields:
                # Search only in specific checkboxes provided by user
                # Map short names (e.g., 'Compound') to full column names ('fields.Compound')
                search_cols = [f"fields.{f}" for f in params.selected_fields if f"fields.{f}" in DF.columns]
                if search_cols:
                    # Check if query string exists in ANY of the selected columns row-wise
                    mask &= DF[search_cols].apply(lambda x: x.astype(str).str.lower().str.contains(params.query)).any(axis=1)
            else:
                # Default Search: Compound OR Smiles
                c_match = DF[COL_MAP['compound']].str.lower().str.contains(params.query, na=False)
                s_match = DF[COL_MAP['smiles']].str.lower().str.contains(params.query, na=False)
                mask &= (c_match | s_match)

        # --- D. Pagination ---
        filtered_df = DF[mask]
        total_items = len(filtered_df)
        total_pages = (total_items + params.per_page - 1) // params.per_page
        
        # Proper slicing
        start = (params.page - 1) * params.per_page
        end = start + params.per_page
        paginated_df = filtered_df.iloc[start:end]

        # --- E. Reconstruct JSON Response ---
        # Convert DataFrame rows back to the nested dictionary format expected by frontend
        items = []
        for _, row in paginated_df.iterrows():
            # Extract keys starting with 'fields.' back into a 'fields' dict
            fields_data = {k.replace('fields.', ''): v for k, v in row.items() if k.startswith('fields.') and pd.notna(v)}
            items.append({
                "mol_idx": int(row['mol_idx']),
                "fields": fields_data,
                "structure_image": row.get('structure_image', None),
                "base64_png": row.get('base64_png', None)
            })

        return jsonify({
            "items": items,
            "page": params.page,
            "per_page": params.per_page,
            "total_pages": total_pages,
            "total_items": total_items
        })

    except Exception as e:
        logger.error("Error in /cytochal-api/items: %s", str(e), exc_info=True)
        return jsonify({"error": "Internal server error"}), 500

@app.get(f"/{CYTOCHAL_API_BASE_PATH}/item/<int:mol_idx>")
def api_item_detail(mol_idx: int):
    try:
        # Find the row where mol_idx matches
        # index[0] gives us the integer position (0, 1, 2...) in the DataFrame
        matches = DF.index[DF['mol_idx'] == mol_idx].tolist()
        
        if not matches:
            return jsonify({"error": "Compound not found"}), 404
            
        current_loc = matches[0]
        
        # Extract Item Data
        row = DF.iloc[current_loc]
        fields_data = {k.replace('fields.', ''): v for k, v in row.items() if k.startswith('fields.') and pd.notna(v)}
        
        # Calculate Next/Prev using iloc (integer location)
        # We rely on the fact that DF was sorted by mol_idx at startup
        prev_idx = int(DF.iloc[current_loc - 1]['mol_idx']) if current_loc > 0 else None
        next_idx = int(DF.iloc[current_loc + 1]['mol_idx']) if current_loc < len(DF) - 1 else None

        return jsonify({
            "item": {
                "mol_idx": int(row['mol_idx']),
                "fields": fields_data,
                "base64_png": row.get('base64_png') or row.get('structure_image'),
            },
            "prev_idx": prev_idx,
            "next_idx": next_idx,
        })

    except Exception as e:
        logger.error("Error in /cytochal-api/item/%s: %s", mol_idx, str(e), exc_info=True)
        return jsonify({"error": "Internal server error"}), 500

if __name__ == "__main__":
    app.run(
        host='0.0.0.0', 
        port=BACKEND_PORT, 
        debug=DEBUG
    )