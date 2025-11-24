# params.py (updated)
from dataclasses import dataclass
from typing import Optional, List

@dataclass
class FilterParams:
    """Dataclass to hold all filter parameters from frontend"""
    query: str = ""
    page: int = 1
    per_page: Optional[int] = None
    min_molweight: Optional[float] = None
    max_molweight: Optional[float] = None
    min_ic50: Optional[float] = None
    max_ic50: Optional[float] = None
    activity: Optional[str] = None
    reversibility: Optional[str] = None
    quantity_type: Optional[str] = None
    min_quantity: Optional[float] = None
    max_quantity: Optional[float] = None
    selected_fields: List[str] = None

    def __post_init__(self):
        if self.selected_fields is None:
            self.selected_fields = []
        if not isinstance(self.selected_fields, list):
            self.selected_fields = []