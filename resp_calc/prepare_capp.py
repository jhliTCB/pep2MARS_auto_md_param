#!/public/home/lijunhao/soft/conda_envs/p4env/bin/python
#TODO: when the covalent residue is otehrs, like SER,LYS, etc...
from pymol import cmd
import sys
import os

covResDict = {'CYS':'CYI','GLU':'GUU','SER':'SRR','LYS':'LYY',
              'HIS':'HII','HIP':'HII','HIE':'HII','HID':'HII',
              'ASP':'ASS','TYR':'TYY','THR':'THH'}

nonLigs = ["NMA","NME","ACE","1PE","2HT","2PE","7PE","ACT","ACY","AKG","BCT",
           "BMA","BME","BOG","BU3","BUD","CAC","CIT","CME","CO3","DMS","DTT",
           "DTV","EDO","EGL","EPE","FES","FMT","FS3","FS4","GBL","GOL","GSH",
           "HEC","HED","HEM","IMD","IOD","IPA","MAN","MES","MG8","MLI","MO6",
           "MPD","MYR","NAG","NCO","NH3","NO3","OCT","OGA","OPG","P2U","PG4",
           "PGE","PGO","PHO","PLP","PO4","POP","PSE","PSU","PTL","SEO","SGM",
           "SO4","SPD","SPM","SRT","SUC","SUL","TAM","TAR","TFA","TLA","TPP",
           "PEG","EOH","BTB","TRS","CAE"]

#TODO: Currently, only single chain is supported, further supports for homo- and
# hetero- oligomers are needed!
def prepare_capp(pdb,ligSele):
    #cmd.reinitialize()
    cmd.load(pdb, "ligress")
    cmd.load(pdb, "cov_md_comp")
    if ligSele == "organic":
        ligSele_string = f"organic & ! r. {'+'.join(nonLigs)}"
    else:
        ligSele_string = ligSele

    cmd.select("cov_organic", f"{ligSele_string} & cov_md_comp")
    ligResName = cmd.get_model("cov_organic").atom[0].resn
    #cmd.select("preCovRES", "cov_organic a. 2.0 & polymer.protein & cov_md_comp")
    cmd.select("preCovRES", "(cov_organic & ! elem H) a. 2.0 & polymer.protein & cov_md_comp")
    covResList = []
    if len(cmd.get_model("preCovRES").atom) > 0:
        for x in cmd.get_model("preCovRES").atom:
            if not x.name.startswith("H"):
                covResList.append((x.resn,x.resi,x.name))
    if len(covResList) <= 1:
        covRes = covResList[0]
    else:
        print('Wanring! more than one residues are found to have non-hydrogen atoms being close\n\
        (<2.0 angstrom) to the ligand, I will exclude those pairs with hydrogen and \n\
        try find the smallest one that is > 1.0 angstrom')
        dists = []
        for i in range(0,len(covResList)):
            dist = min([cmd.get_distance(covResList[i][2],x.name)
                        for x in cmd.get_model("cov_organic").atom])
            dists.append(dist)
            if i > 0 and 1.0 < dist < min(dists[0:i]):
                covRes = covResList[i]
    # covRes [resn, resi, name]
    covResName = covRes[0]
    cmd.select("CovRES", f"r. {covResName} & i. {covRes[1]} & cov_md_comp")
    cmd.alter("CovRES", f"resn='{covResDict[covResName]}'")
    cmd.select("cov_md_comp_organic", f"{ligSele_string} & cov_md_comp" )
    cmd.remove("elem H & cov_md_comp")
    if not cmd.get_model("cov_md_comp_organic").atom[0].chain:
        cmd.alter("cov_md_comp_organic", "chain='X'")
    cmd.set("retain_order", 1)
    cmd.save("covComPdb4tleap.pdb", "cov_md_comp")
    cmd.select("CovRES2", f"r. {covResName} & i. {covRes[1]} & ligress")
    cmd.select("cov_organic_2", f"r. {ligResName} & ligress")
    cmd.select("ace", "(byres (CovRES2 & name N) a. 1.5 & ! CovRES2) & name CA+C+O & ligress")
    cmd.select("nme", "(byres (CovRES2 & name C) a. 1.5 & ! CovRES2) & name N+H+CA & ligress")
    cmd.select("SG", "r. {} & i. {} & name {} & ligress".format(*covRes))
    cmd.select("beta", "SG a. 2.05 & cov_organic_2 & ligress")
    beta = cmd.get_model("beta").atom[0]
    bondinfo = f"REMARK bondinfo: p.LIS.{covRes[2]} p.LIS.{beta.name}\n"
    leapping = f"REMARK bondinfo: cpx.{covResDict[covResName]}.{covRes[2]} cpx.{ligResName}.{beta.name}"
    os.system(f'''sed -i "1i{leapping}" covComPdb4tleap.pdb''')
    cmd.remove("ligress & !(cov_organic_2 | CovRES2 | nme | ace)")
    cmd.alter("ace", "resn='ACE'")
    cmd.alter("nme", "resn='NME'")
    cmd.h_add(selection="r. ACE+NME & ligress")
    cmd.select("mod_h", "(name CA & r. ACE+NME) a. 1.5 & elem H")
    cmd.alter("mod_h","name='HA'")
    #cmd.save("capped_lig_for_esp.pdb")
    cmd.save("capped_residue.pdb", "(r. NME+ACE|CovRES2) & ligress")
    #cmd.set("retain_order",1)
    cmd.save("lig_retained.pdb", "cov_organic_2 & !(r. NME+ACE) & ligress")
    with open("automated_lig_capped.pdb",'w') as fp:
        fp.write(bondinfo)
        for line in open("capped_residue.pdb", "r"):
            if line[0:6] in ['ATOM  ','HETATM']:
                fp.write(line)
        for line in open("lig_retained.pdb",'r'):
            if line[0:6] in ['ATOM  ','HETATM']:
                fp.write(line)
        fp.write('END\n')
    cmd.reinitialize()
    return ["covComPdb4tleap.pdb","capped_residue.pdb","lig_retained.pdb","automated_lig_capped.pdb"]

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"{sys.argv[0]} prepared_covPose.pdb ligsele")
        exit(1)
    elif len(sys.argv) == 2:
        pdbList = prepare_capp(pdb=sys.argv[1],ligSele="organic")
    elif len(sys.argv) == 3:
        pdbList = prepare_capp(pdb=sys.argv[1],ligSele=sys.argv[2])

    for pdb in pdbList:
        print(f"Renaming CL to Cl for amb for {pdb}!")
        os.system(f"sed -i 's/CL/Cl/g' {pdb}")
