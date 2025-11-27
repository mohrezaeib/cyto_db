# parse_sdf.py
import base64
import io
import json
import unicodedata
import difflib
import os
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor
from dotenv import load_dotenv


#load the .env 
load_dotenv()


rdDepictor.SetPreferCoordGen(True)

# Get paths from environment variables with fallbacks
BASE_DIR = Path(__file__).parent.parent
SDF_PATH = BASE_DIR / os.getenv('CYTOCHAL_SDF_PATH', 'static/CytoLabs_Database.sdf')
IMAGE_DIR = BASE_DIR / os.getenv('CYTOCHAL_IMAGE_DIR', 'static/microscopy_images')
STRUCTURE_IMAGE_DIR = BASE_DIR / os.getenv('CYTOCHAL_STRUCTURE_DIR', 'static/molecule_structures_image')

IMAGE_EXTS = {".tif", ".tiff", ".png", ".jpg", ".jpeg"}

def normalize(s: str) -> str:
    """Lowercase, strip accents, and remove non-alphanumerics for robust matching."""
    if not isinstance(s, str):
        s = str(s)
    s = unicodedata.normalize("NFKD", s).casefold()
    return "".join(ch for ch in s if ch.isalnum())

def index_images(image_dir=IMAGE_DIR):
    """Return [(filename, normalized_stem)] for candidates in the image directory."""
    image_dir = Path(image_dir)
    candidates = []
    if image_dir.is_dir():
        for p in image_dir.iterdir():
            if p.is_file() and p.suffix.lower() in IMAGE_EXTS:
                candidates.append((p.name, normalize(p.stem)))
    else:
        print(f"WARNING: Image directory {image_dir} does not exist")
    return candidates

def match_image(compound: str, candidates, threshold: float = 0.80):
    """Find best filename for the given compound, or None if no good match."""
    if not compound:
        return None
    c_norm = normalize(compound)
    best_score = 0.0
    best_name = None
    for fname, norm in candidates:
        if c_norm == norm:
            score = 1.0
        elif c_norm in norm or norm in c_norm:
            score = 0.95
        else:
            score = difflib.SequenceMatcher(None, c_norm, norm).ratio()
        if score > best_score:
            best_score, best_name = score, fname
    return best_name if best_score >= threshold else None

from PIL import Image, ImageOps
def ensure_web_image(path: Path) -> Path:
    if path.suffix.lower() in {".tif", ".tiff"}:
        out = path.with_suffix(".png")
        if (not out.exists()) or (path.stat().st_mtime > out.stat().st_mtime):
            with Image.open(path) as im:
                if getattr(im, "n_frames", 1) > 1:
                    im.seek(0)
                im = ImageOps.exif_transpose(im)
                if im.mode not in ("RGB", "RGBA"):
                    im = im.convert("RGBA")
                out.parent.mkdir(parents=True, exist_ok=True)
                im.save(out, "PNG", optimize=True)
        return out
    return path

def convert_string_to_number(value):
    if value is None or value.strip() == "":
        return ""
    
    value = str(value).strip()
    
    if any(char in value for char in ['/', '\\', '@', '(', ')', '[', ']', '=', '#', '+', '-']):
        return value
    
    import re
    quantity_pattern = r'^\s*[-+]?\d*\.?\d+\s*(mg|g|ml|l|mM|µM|nM|ppm|%|units?)\s*$'
    if re.match(quantity_pattern, value, re.IGNORECASE):
        number_match = re.search(r'[-+]?\d*\.?\d+', value)
        if number_match:
            try:
                num_str = number_match.group()
                f = float(num_str)
                if f.is_integer():
                    return int(f)
                return f
            except ValueError:
                return value
    
    try:
        f = float(value)
        if f.is_integer():
            return int(f)
        return f
    except ValueError:
        return value

def parse_sdf_to_data(sdf_path=SDF_PATH, image_dir=IMAGE_DIR):
    sdf_path = Path(sdf_path)
    if not sdf_path.exists():
        print(f"ERROR: SDF file not found at {sdf_path}")
        print(f"Current working directory: {os.getcwd()}")
        return []
    
    candidates = index_images(image_dir)
    suppl = Chem.SDMolSupplier(str(sdf_path), sanitize=False, removeHs=False, strictParsing=False)
    
    STRUCTURE_IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    
    first_mol = next(iter(suppl), None)
    if first_mol is not None:
        print("First molecule properties:", list(first_mol.GetPropNames()))
    else:
        print("No valid molecules found in SDF")

    results = []
    mol_idx = 0

    for mol in suppl:
        if mol is None:
            print("Skipping invalid molecule")
            continue

        Chem.rdDepictor.Compute2DCoords(mol)
        img = Draw.MolToImage(mol, size=(1000, 1000))

        STRUCTURE_IMAGE_DIR.mkdir(parents=True, exist_ok=True)

        structure_filename = f"structure_{mol_idx}.png"
        structure_path = STRUCTURE_IMAGE_DIR / structure_filename
        img.save(structure_path, format='PNG')

        structure_file_ref = structure_filename  

        fields_dict = {}
        for prop_name in mol.GetPropNames():
            value = mol.GetProp(prop_name)
            converted = convert_string_to_number(value)

            if prop_name == "Reversibilty" and converted == "":
                converted = "not tested"
            elif prop_name == "Actin Disruption Activity" and converted == "":
                converted = "not tested"
            elif prop_name == "Quantity" and converted == "":
                converted = "not available"

            fields_dict[prop_name] = converted

        compound_name = fields_dict.get("Compound", "")
        matched = match_image(compound_name, candidates)
        if matched:
            p = (IMAGE_DIR / matched).resolve()
            p = ensure_web_image(p)
            fields_dict["Image File"] = p.name

        results.append({
            "mol_idx": mol_idx,
            "structure_image": structure_file_ref,
            "fields": fields_dict,
        })
        mol_idx += 1

    return results

if __name__ == "__main__":
    print(f"Using SDF file: {SDF_PATH}")
    print(f"SDF file exists: {SDF_PATH.exists()}")
    
    data = parse_sdf_to_data()
    
    output_path = BASE_DIR / "static" / "data.json"
    with open(output_path, "w", encoding="utf-8") as file:
        json.dump(data, file, ensure_ascii=False, indent=2)
    print(f"Data written to: {output_path}")