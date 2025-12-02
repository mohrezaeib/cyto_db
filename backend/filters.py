# filters.py (updated with pandas)
import re
from typing import Any, Dict, List, Optional
import pandas as pd
from search_engine import search_engine
from params import FilterParams

class CompoundFilter:
    """Main filter class for applying all filters to compound data using pandas"""

    def __init__(self, cache=None):
        self._num_re = re.compile(r"[-+]?\d*\.\d+|\d+")
        self.cache = cache

    def get_field_value(self, row: pd.Series, target_field: str) -> Any:
        """Get field value with flexible naming matching from pandas Series"""
        target_lower = target_field.lower().replace(" ", "")
        
        # Check if target_field exists directly in row index
        if target_field in row.index:
            return row[target_field]
        
        # Check fields dictionary if it exists
        fields_dict = row.get('fields', {})
        if isinstance(fields_dict, dict):
            for field_name, value in fields_dict.items():
                field_lower = field_name.lower().replace(" ", "")
                if target_lower == field_lower:
                    return value
                if target_lower in field_lower or field_lower in target_lower:
                    return value
        
        # Check direct row fields
        for field_name in row.index:
            field_lower = str(field_name).lower().replace(" ", "")
            if target_lower == field_lower:
                return row[field_name]
        
        return None

    def parse_number(self, value: Any) -> Optional[float]:
        """Extract first numeric value from any input"""
        if value is None:
            return None
        if isinstance(value, (int, float)):
            return float(value)
        if pd.isna(value):
            return None
        match = self._num_re.search(str(value))
        return float(match.group()) if match else None

    def apply_numeric_filter(self, item_value: Any, min_val: Optional[float], max_val: Optional[float]) -> bool:
        """Apply min/max numeric filter with None handling"""
        numeric_value = self.parse_number(item_value)
        if min_val is None and max_val is None:
            return True
        if numeric_value is None:
            return False
        if min_val is not None and numeric_value < min_val:
            return False
        if max_val is not None and numeric_value > max_val:
            return False
        return True

    def apply_activity_filter(self, row: pd.Series, activity: str) -> bool:
        """Apply actin disruption activity filter in a simple, robust way."""
        if not activity:
            return True

        fields = row.get("fields", {}) or {}

        activity_value = ""
        for name, value in fields.items():
            key_norm = name.lower()
            if "actin disruption activity".lower() in key_norm:
                activity_value = str(value or "")
                break

        if not activity_value:
            # no matching field found for this row
            return False

        # e.g. activity="++" → matches "++", "+/-", etc. as you like
        return activity.lower() in activity_value.lower()


    def apply_reversibility_filter(self, row: pd.Series, reversibility: str) -> bool:
        """Apply reversibility filter with mapping"""
        if not reversibility:
            return True
        
        reversibility_value = str(
            self.get_field_value(row, "Reversibilty") or 
            self.get_field_value(row, "Reversibility") or ""
        )
        
        if reversibility == "+":
            return "+" in reversibility_value
        elif reversibility == "-":
            return "-" in reversibility_value
        elif reversibility.lower() == "not tested":
            return "not tested" in reversibility_value.lower()
        else:
            return reversibility.lower() in reversibility_value.lower()

    def apply_quantity_filter(self, row: pd.Series, params: FilterParams) -> bool:
        """Apply quantity filter based on quantity_type"""
        if not params.quantity_type:
            return True
        
        quantity_value = self.get_field_value(row, "quantity")
        quantity_numeric = self.parse_number(quantity_value)
        
        if params.quantity_type == "numeric":
            return self.apply_numeric_filter(quantity_value, params.min_quantity, params.max_quantity)
        elif params.quantity_type == "available":
            return str(quantity_value).lower() == "available"
        elif params.quantity_type == "not available":
            return not (str(quantity_value).lower() == "available" or quantity_numeric is not None)
        return True

    def apply_search_filter(self, row: pd.Series, params: FilterParams) -> bool:
        """Apply text search using Whoosh with Flask caching"""
        if not params.query or not params.query.strip():
            return True
        
        # Create cache key from params
        cache_key = f"search_{params.query}_{'_'.join(sorted(params.selected_fields))}"
        
        # Try to get from cache first
        if self.cache:
            cached_result = self.cache.get(cache_key)
            if cached_result is not None:
                doc_id = str(row.get("mol_idx") or hash(str(row.to_dict())))
                return doc_id in cached_result
        
        # Get document ID from the row
        doc_id = str(row.get("mol_idx") or hash(str(row.to_dict())))
        
        # Use exact field names that match your JSON
        field_map = {
            "Compound": "Compound",
            "Origin": "Origin", 
            "Reference": "Reference",
            "Smiles": "Smiles", 
            "Subclass": "Subclass",
            "Molecular Formula": "Molecular_Formula"
        }
        
        # Convert selected fields to Whoosh fields
        whoosh_fields = []
        if params.selected_fields:
            whoosh_fields = [field_map[f] for f in params.selected_fields if f in field_map]
        else:
            whoosh_fields = list(field_map.values())
        
        # Search and check if this document matches
        matching_ids = search_engine.search(params.query.strip(), whoosh_fields)
        
        # Cache the result
        if self.cache:
            self.cache.set(cache_key, matching_ids, timeout=300)
        
        return doc_id in matching_ids

    def filter_compounds(self, compounds_df: pd.DataFrame, params: FilterParams) -> pd.DataFrame:
        """Main filtering method with pandas"""
        if compounds_df.empty:
            return compounds_df
        
        # Create cache key from all parameters
        cache_key = f"filter_{params.query}_{params.page}_{params.per_page}_{params.min_molweight}_{params.max_molweight}_{params.min_ic50}_{params.max_ic50}_{params.activity}_{params.reversibility}_{params.quantity_type}_{params.min_quantity}_{params.max_quantity}_{'_'.join(sorted(params.selected_fields))}"
        
        # Try to get from cache first
        if self.cache:
            cached_result = self.cache.get(cache_key)
            if cached_result is not None:
                return cached_result
        
        # Start with all compounds
        filtered_df = compounds_df.copy()
        
        # Apply filters using pandas boolean indexing
        mask = pd.Series([True] * len(filtered_df), index=filtered_df.index)
        
        # Apply numeric filters
        if params.min_molweight is not None or params.max_molweight is not None:
            molweight_mask = filtered_df.apply(
                lambda row: self.apply_numeric_filter(
                    self.get_field_value(row, "totalmolweight"),
                    params.min_molweight,
                    params.max_molweight
                ), axis=1
            )
            mask = mask & molweight_mask
        
        if params.min_ic50 is not None or params.max_ic50 is not None:
            ic50_mask = filtered_df.apply(
                lambda row: self.apply_numeric_filter(
                    self.get_field_value(row, "ic50"),
                    params.min_ic50,
                    params.max_ic50
                ), axis=1
            )
            mask = mask & ic50_mask
        
        # Apply activity filter
        if params.activity:
            activity_mask = filtered_df.apply(
                lambda row: self.apply_activity_filter(row, params.activity), axis=1
            )
            mask = mask & activity_mask
        
        # Apply reversibility filter
        if params.reversibility:
            reversibility_mask = filtered_df.apply(
                lambda row: self.apply_reversibility_filter(row, params.reversibility), axis=1
            )
            mask = mask & reversibility_mask
        
        # Apply quantity filter
        if params.quantity_type:
            quantity_mask = filtered_df.apply(
                lambda row: self.apply_quantity_filter(row, params), axis=1
            )
            mask = mask & quantity_mask
        
        # Apply search filter
        if params.query and params.query.strip():
            search_mask = filtered_df.apply(
                lambda row: self.apply_search_filter(row, params), axis=1
            )
            mask = mask & search_mask
        
        # Apply the final mask
        filtered_df = filtered_df[mask]
        
        # Cache the filtered result
        if self.cache:
            self.cache.set(cache_key, filtered_df, timeout=300)
        
        return filtered_df