#!/public/home/lijunhao/soft/conda_envs/p4env/bin/python
import sys, os
from pymol import cmd
from pymol import stored
import subprocess

def add_numbers_to_elements(lst):
    counts = {}
    result = []
    for element in lst:
        counts[element] = counts.get(element, 0) + 1
        result.append(f"{element}{counts[element]}")
    return result

def preparePDB(cappedPdb,ligCharge):
    cmd.load(cappedPdb)
    cmd.set("retain_order",1)
    resNames = {"ACE":[x.name for x in cmd.get_model("r. ACE").atom],
                "NME":[x.name for x in cmd.get_model("r. NME").atom],
                "CYS":[x.name for x in cmd.get_model("r. CYS").atom]}
    # unify the atom names
    stored.elements = []
    cmd.do(f"iterate all, stored.elements.append(elem)")
    #elements = [x.elem.upper() for x in cmd.get_model("all").atom]
    newNames = add_numbers_to_elements(stored.elements) #unified names
    allAtoms = cmd.get_model("all").atom
    oldNames = [allAtoms[x].name for x in range(len(allAtoms))]
    print(f"checking: {newNames}")
    i = 0
    for atom in allAtoms:
        cmd.alter(f"rank {i}",f"name='{newNames[i]}'")
        i += 1

    #check namings:
    altered_names = [x.name for x in cmd.get_model("all").atom]
    if altered_names != newNames:
        print("something wrong!")
        print(f"newNames {newNames}")
        print(f"alteredNames {altered_names}")
    #update names:
    uniNames = {"ACE": [x.name for x in cmd.get_model("r. ACE").atom],
                "NME": [x.name for x in cmd.get_model("r. NME").atom],
                "CYS":[x.name for x in cmd.get_model("r. CYS").atom]}
    omitName = set(uniNames["ACE"]+uniNames["NME"])
    #write the omiting atom list for prepgen
    with open("lig.mc",'w') as fp:
        fp.write(f"HEAD_NAME {uniNames['CYS'][resNames['CYS'].index('N')]}\n")
        fp.write(f"TAIL_NAME {uniNames['CYS'][resNames['CYS'].index('C')]}\n")
        fp.write(f"MAIN_CHAIN {uniNames['CYS'][resNames['CYS'].index('CA')]}\n")
        for omi in omitName:
            fp.write(f"OMIT_NAME {omi}\n")
        fp.write(f"PRE_HEAD_TYPE C\n")
        fp.write(f"POST_TAIL_TYPE N\n")
        fp.write(f"CHARGE {ligCharge}\n")
    i = 1
    #unify others
    for atom in cmd.get_model("all").atom:
        cmd.alter(f"name {atom.name}","resn='LIS'")
        cmd.alter(f"name {atom.name}","resi='999'")
        cmd.alter(f"name {atom.name}","chain='A'")
        cmd.alter(f"name {atom.name}",f"id='{i}'")
        i += 1
    cmd.set("retain_order", 0)
    cmd.save("cappedPdb.pdb")
    sed_cmd = ["sed","-i","""'/CONECT/d;s/ATOM  /HETATM/g;/TER/d'""","cappedPdb.pdb"]
    print(f"sed_cmd is {' '.join(sed_cmd)}")
    os.system(' '.join(sed_cmd))
    with open("newAtomNames.lst","w") as fp: fp.write('\n'.join(newNames))
    with open("oldAtomNames.lst","w") as fp: fp.write('\n'.join(oldNames))
    return ["cappedPdb.pdb", oldNames, newNames]

def fit_covONE_bcc(pdb,ligCharge):
    _cmd1 = ["antechamber","-i",pdb,"-fi","pdb","-o","lis.ac",
            "-c","bcc","-s","2","-at","amber","-nc",str(ligCharge),
            "-dr","n","-rn","LIS","-fo","ac"]
    print(f"running: {' '.join(_cmd1)}")
    bcc = subprocess.run(_cmd1,capture_output=True,text=True)
    _cmd2 = ["prepgen","-i","lis.ac","-f","int",
             "-m","lig.mc","-rn","LIS","-o","lis.prepi"]
    print(f"running: {' '.join(_cmd2)}")
    pre = subprocess.run(_cmd2,capture_output=True,text=True)
    _cmd3 = ["parmchk2","-i","lis.prepi","-f","prepi","-o","lis.frcmod","-a","Y"]
    print(f"running: {' '.join(_cmd3)}")
    frc = subprocess.run(_cmd3,capture_output=True,text=True)
    return bcc+pre+frc

def fit_covONE_resp(pdbFile,espFile):
    _cmd1 = ["antechamber", "-i", espFile, "-fi", "gesp", "-o", "lis.ac",
             "-c", "resp", "-s", "2", "-at", "amber",
             "-dr", "n", "-rn", "LIS", "-fo", "ac"]
    print(f"running: {' '.join(_cmd1)}")
    bcc = subprocess.run(_cmd1, capture_output=True, text=True)
    _cmd1 = ["antechamber", "-i", "lis.ac", "-fi", "ac", "-c", "wc", "-cf", "lis.chg"]
    os.system(' '.join(_cmd1))
    _cmd1 = ["antechamber", "-i", pdbFile, "-fi", "pdb", "-o", "lis_pdb.ac",
             "-fo", "ac", "-c", "rc", "-cf", "lis.chg"]
    _cmd2 = ["prepgen", "-i", "lis_pdb.ac", "-f", "int",
             "-m", "lig.mc", "-rn", "LIS", "-o", "lis.prepi"]
    print(f"running: {' '.join(_cmd2)}")
    pre = subprocess.run(_cmd2, capture_output=True, text=True)
    _cmd3 = ["parmchk2", "-i", "lis.prepi", "-f", "prepi", "-o", "lis.frcmod", "-a", "Y"]
    print(f"running: {' '.join(_cmd3)}")
    frc = subprocess.run(_cmd3, capture_output=True, text=True)

if __name__ == "__main__":
    pdbFile = sys.argv[1]
    ligChrg = sys.argv[2]
    espFile = sys.argv[3]
    prepared = preparePDB(pdbFile, ligChrg)[0]
    if espFile == 'bcc':
        results = fit_covONE_bcc(prepared,ligChrg)
    elif espFile.endswith(".esp"):
        results = fit_covONE_resp(pdbFile,espFile)
