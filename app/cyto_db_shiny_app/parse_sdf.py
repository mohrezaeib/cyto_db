# parse_sdf.py
import base64
import io
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor

rdDepictor.SetPreferCoordGen(True)

def parse_sdf_to_data(sdf_path):
    suppl = Chem.SDMolSupplier(sdf_path, sanitize=True, removeHs=True, strictParsing=True)
    results = []
    mol_idx = 0

    for mol in suppl:
        if mol is None:
            continue
        Chem.rdDepictor.Compute2DCoords(mol)
        img = Draw.MolToImage(mol, size=(1000,1000))
        buf = io.BytesIO()
        img.save(buf, format='PNG')
        b64 = base64.b64encode(buf.getvalue()).decode('utf-8')

        fields_dict = {}
        for prop_name in mol.GetPropNames():
            fields_dict[prop_name] = mol.GetProp(prop_name)

        results.append({'mol_idx': mol_idx, 'base64_png': b64, 'fields': fields_dict})
        mol_idx += 1

    return results
