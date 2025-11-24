# filters.py (updated)
import re
from typing import Any, Dict, List, Optional
from search_engine import search_engine
from params import FilterParams

class CompoundFilter:
    """Main filter class for applying all filters to compound data"""

    def __init__(self, cache=None):
        self._num_re = re.compile(r"[-+]?\d*\.\d+|\d+")
        self.cache = cache

    def get_field_value(self, item_fields: Dict[str, Any], target_field: str) -> Any:
        """Get field value with flexible naming matching"""
        target_lower = target_field.lower().replace(" ", "")
        
        for field_name, value in item_fields.items():
            field_lower = field_name.lower().replace(" ", "")
            if target_lower == field_lower:
                return value
            if target_lower in field_lower or field_lower in target_lower:
                return value
        return None

    def parse_number(self, value: Any) -> Optional[float]:
        """Extract first numeric value from any input"""
        if value is None:
            return None
        if isinstance(value, (int, float)):
            return float(value)
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

    def apply_activity_filter(self, fields: Dict[str, Any], activity: str) -> bool:
        """Apply actin disruption activity filter"""
        if not activity:
            return True
        activity_value = str(fields.get("Actin Disruption Activity") or "")
        return activity.lower() in activity_value.lower()

    def apply_reversibility_filter(self, fields: Dict[str, Any], reversibility: str) -> bool:
        """Apply reversibility filter with mapping"""
        if not reversibility:
            return True
        reversibility_value = str(fields.get("Reversibilty") or fields.get("Reversibility") or "")
        if reversibility == "+":
            return "+" in reversibility_value
        elif reversibility == "-":
            return "-" in reversibility_value
        elif reversibility.lower() == "not tested":
            return "not tested" in reversibility_value.lower()
        else:
            return reversibility.lower() in reversibility_value.lower()

    def apply_quantity_filter(self, fields: Dict[str, Any], params: FilterParams) -> bool:
        """Apply quantity filter based on quantity_type"""
        if not params.quantity_type:
            return True
        quantity_value = self.get_field_value(fields, "quantity")
        quantity_numeric = self.parse_number(quantity_value)
        if params.quantity_type == "numeric":
            return self.apply_numeric_filter(quantity_value, params.min_quantity, params.max_quantity)
        elif params.quantity_type == "available":
            return str(quantity_value).lower() == "available"
        elif params.quantity_type == "not available":
            return not (str(quantity_value).lower() == "available" or quantity_numeric is not None)
        return True

    def apply_search_filter(self, compound, params: FilterParams) -> bool:
        """Apply text search using Whoosh with Flask caching"""
        if not params.query or not params.query.strip():
            return True
        
        # Create cache key from params
        cache_key = f"search_{params.query}_{'_'.join(sorted(params.selected_fields))}"
        
        # Try to get from cache first
        if self.cache:
            cached_result = self.cache.get(cache_key)
            if cached_result is not None:
                doc_id = str(compound.get("mol_idx") or hash(str(compound)))
                return doc_id in cached_result
        
        fields = compound.get("fields", {})
        
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
        
        # Get document ID from the compound (not just fields)
        doc_id = str(compound.get("mol_idx") or hash(str(compound)))
        
        # Search and check if this document matches
        matching_ids = search_engine.search(params.query.strip(), whoosh_fields)
        
        # Cache the result
        if self.cache:
            self.cache.set(cache_key, matching_ids, timeout=300)
        
        return doc_id in matching_ids

    def filter_compounds(self, compounds: List[Dict], params: FilterParams) -> List[Dict]:
        """Main filtering method with Flask caching"""
        # Create cache key from all parameters
        cache_key = f"filter_{params.query}_{params.page}_{params.per_page}_{params.min_molweight}_{params.max_molweight}_{params.min_ic50}_{params.max_ic50}_{params.activity}_{params.reversibility}_{params.quantity_type}_{params.min_quantity}_{params.max_quantity}_{'_'.join(sorted(params.selected_fields))}"
        
        # Try to get from cache first
        if self.cache:
            cached_result = self.cache.get(cache_key)
            if cached_result is not None:
                return cached_result
        
        filtered = []
        for compound in compounds:
            fields = compound.get("fields", {})
            if not all([
                self.apply_numeric_filter(
                    self.get_field_value(fields, "totalmolweight"),
                    params.min_molweight,
                    params.max_molweight
                ),
                self.apply_numeric_filter(
                    self.get_field_value(fields, "ic50"),
                    params.min_ic50,
                    params.max_ic50
                ),
                self.apply_activity_filter(fields, params.activity),
                self.apply_reversibility_filter(fields, params.reversibility),
                self.apply_quantity_filter(fields, params),
                self.apply_search_filter(compound, params)
            ]):
                continue
            filtered.append(compound)
        
        # Cache the filtered result
        if self.cache:
            self.cache.set(cache_key, filtered, timeout=300)
        
        return filtered