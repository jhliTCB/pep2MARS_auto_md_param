#!/public/home/lijunhao/soft/conda_envs/p4env/bin/python
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
def smi2sdf(smi, workDir, molname='structureFromSmiles'):
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
    AllChem.UFFOptimizeMolecule(mol)
    outPath = f"{workDir}/{molname}.sdf"
    mol.SetProp("_Name", molname)
    Chem.SDWriter(outPath).write(mol)
    return outPath

smi2sdf(smi=sys.argv[1],workDir=sys.argv[2],molname=sys.argv[3])