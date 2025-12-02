# app.py (updated with pandas)
from flask import Flask, request, jsonify
from flask_cors import CORS
from flask_caching import Cache
import json
import os
import sys 
import logging
import pandas as pd
from params import FilterParams
from filters import CompoundFilter
from search_engine import search_engine
from dotenv import load_dotenv

load_dotenv()

# Load environment variables
BACKEND_PORT = int(os.getenv('CYTOCHAL_BACKEND_PORT', 5000))
CYTOCHAL_API_BASE_PATH = os.getenv("CYTOCHAL_API_BASE_PATH", "cytochal-api")
FLASK_ENV = os.getenv('CYTOCHAL_FLASK_ENV', 'production')
DEBUG = os.getenv('CYTOCHAL_DEBUG', 'False').lower() == 'true'
DATA_PATH = os.getenv('CYTOCHAL_DATA_PATH', 'static/data.json')
LOG_LEVEL = os.getenv('CYTOCHAL_LOG_LEVEL', 'INFO')

# Configure logging
log_level = getattr(logging, LOG_LEVEL.upper(), logging.INFO)
logging.basicConfig(level=log_level)
logger = logging.getLogger(__name__)

app = Flask(__name__)
CORS(app)

# Configure Flask-Caching
app.config['CACHE_TYPE'] = 'SimpleCache'
app.config['CACHE_DEFAULT_TIMEOUT'] = 300  # 5 minutes
cache = Cache(app)

# Initialize filter with Flask cache
compound_filter = CompoundFilter(cache=cache)

# Global DataFrame to store compounds
compounds_df = None

try:
    if not os.path.exists(DATA_PATH):
        raise FileNotFoundError(f"Data file not found at: {DATA_PATH}")
    
    # If file exists but is truly empty (0 bytes)
    if os.path.getsize(DATA_PATH) == 0:
        raise ValueError(f"Data file at {DATA_PATH} is empty (0 bytes).")

    with open(DATA_PATH, "r", encoding="utf-8") as f:
        data = json.load(f)
        
    # If JSON is valid but contains no records: [] or {}
    if not data:
        raise ValueError(f"Data file at {DATA_PATH} contains no records.")
    
    # Convert to pandas DataFrame
    compounds_df = pd.DataFrame(data)
    
    # If DataFrame ended up empty for any reason
    if compounds_df.empty:
        raise ValueError(f"Data file at {DATA_PATH} produced an empty DataFrame.")
    
    # Ensure mol_idx is set properly
    if 'mol_idx' not in compounds_df.columns:
        compounds_df['mol_idx'] = compounds_df.index
    
    # Index compounds for search
    for _, compound in compounds_df.iterrows():
        doc_id = str(compound.get("mol_idx") or hash(str(compound.get("fields", {}))))
        search_engine.index_document(doc_id, compound.get("fields", {}))
        
    logger.info(f"Successfully loaded {len(compounds_df)} compounds from {DATA_PATH}")
        
except Exception as e:
    # Log full traceback and abort startup
    logger.exception(f"Error loading data: {e}")
    sys.exit(1)

@app.route("/cytochal-api/health")
def health_check():
    return jsonify({
        "status": "healthy", 
        "data_loaded": len(compounds_df),
        "environment": FLASK_ENV
    })

@app.route(f"/{CYTOCHAL_API_BASE_PATH}/items")
def api_items():
    try:
        # Parse query parameters
        params = FilterParams(
            query=request.args.get("query", "").strip(),
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

        # Apply filters using pandas
        filtered_data = compound_filter.filter_compounds(compounds_df, params)

        # Convert filtered DataFrame to list of dictionaries
        filtered_list = filtered_data.to_dict('records')
        
        # Paginate results
        total_items = len(filtered_list)
        total_pages = (total_items + params.per_page - 1) // params.per_page
        start_index = (params.page - 1) * params.per_page
        end_index = start_index + params.per_page
        paginated_items = filtered_list[start_index:end_index]

        result = {
            "items": paginated_items,
            "page": params.page,
            "per_page": params.per_page,
            "total_pages": total_pages,
            "total_items": total_items
        }
        
        return jsonify(result)

    except Exception as e:
        logger.error("Error in /cytochal-api/items: %s", str(e))
        return jsonify({"error": "Internal server error"}), 500

@app.get(f"/{CYTOCHAL_API_BASE_PATH}/item/<int:mol_idx>")
def api_item_detail(mol_idx: int):
    try:
        # Find item using pandas
        item_row = compounds_df[compounds_df['mol_idx'] == mol_idx]
        if item_row.empty:
            return jsonify({"error": "Compound not found"}), 404

        item = item_row.iloc[0].to_dict()
        
        # Find previous/next indexes
        current_index = compounds_df[compounds_df['mol_idx'] == mol_idx].index[0]

        prev_idx = compounds_df.iloc[current_index - 1]["mol_idx"] if current_index > 0 else None
        next_idx = compounds_df.iloc[current_index + 1]["mol_idx"] if current_index < len(compounds_df) - 1 else None

        return jsonify({
            "item": {
                "mol_idx": item.get("mol_idx"),
                "fields": item.get("fields", {}),
                "base64_png": item.get("base64_png"),
            },
            "prev_idx": prev_idx,
            "next_idx": next_idx,
        })

    except Exception as e:
        logger.error("Error in /cytochal-api/item/%s: %s", mol_idx, str(e))
        return jsonify({"error": "Internal server error"}), 500

# Cache management endpoints
@app.route("/cytochal-api/cache/clear", methods=['POST'])
def clear_cache():
    """Clear the entire cache"""
    cache.clear()
    return jsonify({"status": "cache cleared"})

@app.route("/cytochal-api/cache/stats")
def cache_stats():
    """Get cache statistics (simplified for Flask-Caching)"""
    # Flask SimpleCache doesn't provide detailed stats
    return jsonify({"message": "Stats not available for SimpleCache"})

if __name__ == "__main__":
    app.run(
        host='0.0.0.0', 
        port=BACKEND_PORT, 
        debug=DEBUG
    )