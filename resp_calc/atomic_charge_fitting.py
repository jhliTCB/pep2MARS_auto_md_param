#!/public/home/lijunhao/soft/conda_envs/p4dev/bin/python
#https://github.com/cdsgroup/resp/tree/master/examples
# ESP/RESP/BCC charge fitting and Amber force field file generation
# resp package need to be patched (in the p4env)
# by Junhao Li

import argparse
import time
import psi4
import resp
import os
import subprocess
import logging
from pymol import cmd
###from rdkit.Contrib.mmpa.create_mmp_db import results

from pdb2com4resp import gen_gau_inp
from prepare_capp import prepare_capp
#from scripts.toolkits.tool1_pdb_manipulation.pdb_tofasta import pdb_to_fasta

## Global settings
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

exec(open('/public/software/apps/module/init/python.py').read())
module('purge')
module('load amber/a22t23 gaussian/16.c02')
module('load mpi/openmpi/gnu/4.0.3 compiler/cuda/11.8')
module('list')

def read_par():
    # Set default values
    scrp_dir = os.path.dirname(os.path.realpath(__file__))
    pdb_type = 'lig'
    lig_sele = 'organic'
    chg_type = 'resp'
    opt_type = 'none'
    esp_type = 'psi4'
    esp_func = 'hf'
    esp_basis = 'aug-cc-pvtz'
    parser = argparse.ArgumentParser(
        description='Script for atomic charge fitting', 
        formatter_class=argparse.RawTextHelpFormatter,
        usage='%(prog)s [options]'
    )    
    # Add arguments with formatted help text
    parser.add_argument('-wd', '--wdir', dest='workDir',
        help='Working directory for calculations    (Required)\n')
        
    parser.add_argument('-chg', '--lig_charge', dest='lig_charge',
        help='Charge of ligand (non-covalent form)    (Required)'
             f'e.g. -chg 1 or -chg -7)\n')
        
    parser.add_argument('-inp', '--inp_file', dest='inp_file',
        help='Input file (.pdb/.xyz/.sdf) or SMILES string    (Required)\n')
        
    parser.add_argument('-pt', '--pdb_type', dest='pdb_type', default=pdb_type,
        help=f'PDB file type, if input ends with .pdb    Available options:\n'
             f'cov : covalent docked complex\n'
             f'lig : pure ligand file\n'
             f'com : non-covalent ligand/receptor complex\n'
             f'(Default: {pdb_type}, when input ends with .pdb)\n')
        
    parser.add_argument('-ligSel', dest='ligSel', default=lig_sele,
        help=f'Ligand selection text (if -pt is com or cov)    (Default: {lig_sele})\n'
             f'Example: -ligSel B:IAQ:442 (chainName:resName:resID)\n')
        
    parser.add_argument('-ct', '--charge_type', dest='charge_type', default=chg_type,
        help=f'Charge fitting method: resp/bcc/esp\n  (Default: {chg_type})\n')
        
    parser.add_argument('-opt', '--optimization', dest='optimization', default=opt_type,
        help=f'Geometry optimization options for resp fitting    Available options:\n'
             f'None|none : do not perform optimization\n'
             f'gau       : optimize using G16 (b3lyp/6-31g**)\n'
             f'psi4      : optimize using psi4 (b3lyp/6-31g**)\n'
             f'(Default: {opt_type}, no geometrical optimizations)\n')
        
    parser.add_argument('--esp_engine', dest='esp_engine', default=esp_type,
        help=f'ESP calculation engine    Available options:\n'
             f'gau|gaussian\n'
             f'psi4\n'
             f'(Default: {esp_type})')
        
    parser.add_argument('--esp_functional', dest='esp_functional', default=esp_func,
        help=f'Functional for ESP calculation  (Default: {esp_func})')
        
    parser.add_argument('--esp_basis', dest='esp_basis', default=esp_basis,
        help=f'Basis set for ESP calculation\n  (Default: {esp_basis})')
    
    parser.add_argument('--resn', dest='resName',type=str,default="LIG",
        help=f'residue name, default is LIG when no resn found in input\n'
             f'else it will the name found in the input file, or from ligSel')

    parser.add_argument('-scrd', '--script_dir', dest='script_dir', default=scrp_dir,
                        help='Script directory location    (Default: location of this script)\n')

    parser.add_argument("--cpus",dest='cpus', default="24",type=str,
                        help='Number of CPUs to use (Default: 24)\n')
    parser.add_argument('--mems',dest='mems', default=96,type=int,
                        help='Memory size in GB (Default: 96)\n')
    parser.add_argument("--just_opt",action="store_true",dest='just_opt',default=False)
    parser.add_argument("--ac_flags",dest="ac_flags",default=" ",type=str,
                        help='''antechamber args, e.g. "@dr no" script will replace @ by -''')
    return parser

def mapping_atom_back(xyz,pdb):
    # mapping the pdb atom name back to TS xyz
    # for PDB fed to antechamber, remove CONECT records
    lines = [x for x in open(xyz,'r')]
    atom_num = int(lines[0].strip('\n'))
    comment = lines[1]
    xyzs = lines[2:]
    i = 0
    news = [f"REMAKR {atom_num} atoms\nREMARK {comment}"]
    for line in open(pdb, 'r'):
        #if line.startswith("TER"): continue
        #if line.startswith('CONECT'): continue
        if line[0:6] not in ['ATOM  ','HETATM']:
            news.append(line)
        elif xyzs[i].split()[0] in line.strip('\n').split()[-1]:
            newxyz = [float(x) for x in xyzs[i].split()[1:]]
            newxyz = "{:8.3f}{:8.3f}{:8.3f}".format(*newxyz)
            oldxyz = line[30:54]
            #if line[17:20].strip() not in ['ACE', 'NME']:
            news.append(line.replace(oldxyz,newxyz))
            i += 1
        else:
            print("FATAL, atoms are not the same, existing...")
            return None
    return news
    #with open(new,'w') as fp:
    #    fp.write(''.join(news))

def smi2sdf(smi, molname='structureFromSmiles'):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
    AllChem.UFFOptimizeMolecule(mol)
    outPath = f"{molname}.sdf"
    mol.SetProp("_Name", molname)
    Chem.SDWriter(outPath).write(mol)
    return outPath

def read_xyz(xyz,chg='0',flag='psi4'):
    string = "{} 1\n".format(chg.replace('e', '-'))
    pxyzs, gxyzs = [], ''
    units = False
    lineN = 0
    for line in open(xyz,'r'):
        lineN += 1
        if lineN > 2:
            string += line
            gxyzs += line
            pxyzs.append([float(x.strip()) for x in line.split()[1:]])
        if "units angstrom" in line:
            units = True
    if not units:
        string += "units angstrom\n"
    if flag == 'psi4': return string
    if flag == 'gau': return gxyzs
    if flag == 'pureXYZ': return pxyzs
   
def psi4_esp_run(inp_file,lig_charge,optimization,functional,basis,cpus,mems):
    basename = inp_file[:-4] #assuming it is a PDB file
    esp_file = f"results_{basename}.esp"
    out = f"{basename}.out"
    psi4.core.set_output_file(out, True)
    psi4_io = psi4.core.IOManager.shared_object()
    psi4_io.set_default_path("/public/cadd/scratch/psi4")
    xyz = f"{basename}.xyz"
    cmd.load(inp_file)
    cmd.set("retain_order",1)
    cmd.save(xyz)
    cmd.reinitialize()
    mol = psi4.geometry(read_xyz(xyz,lig_charge))
    mol.update_geometry()
    mol.set_name(basename)
    psi4.set_memory(f"{mems}GB")
    psi4.set_num_threads(int(cpus))
    if os.path.isfile(esp_file):
        for line in open(esp_file,'r'):
            if "ESP VALUES AND GRID POINT COOR" in line:
                logging.info("esp file exists, skipping psi4 run")
                return esp_file
    if optimization:
        mol.save_xyz_file(f"{basename}_init.xyz",1)
        psi4.optimize("b3lyp/6-31g*", molecule=mol)
        mol.save_xyz_file(f"{basename}_opt.xyz",1)
        cmd.load(f"{basename}_opt.xyz")
        cmd.save(f"{basename}_opt.pdb")
        cmd.reinitialize()
        #ase.io.write(basename+"_opt.pdb", ase.io.read(basename+"_opt.xyz"))
    else:
        mol.save_xyz_file(f"{basename}_init.xyz",1)
        cmd.load(f"{basename}_init.xyz")
        cmd.save(f"{basename}_init.pdb")
        cmd.reinitialize()
        #ase.io.write(basename+"_init.pdb", ase.io.read(basename+"_init.xyz"))
        #write_pdb("just_xyz", xyz[1:], basename+".pdb")

    mol.update_geometry()
    ## calculation of ESP:
    options = {'VDW_SCALE_FACTORS'  : [1.4, 1.6, 1.8, 2.0],
               'VDW_POINT_DENSITY'  : 1.0,
               'RESP_A'             : 0.0005,
               'RESP_B'             : 0.1,
               'METHOD_ESP'         : functional,
               'BASIS_ESP'          : basis,
               'CPU'                : int(cpus),
               'MEM'                : f"{mems}GB"
               }
    charges1 = resp.resp([mol], options)
    return esp_file

def run_gau(basename,cpus):
    #wrapper for gaussian run after gen_gau_inp
    script_dir = os.path.dirname(os.path.realpath(__file__))
    _cmd = [f"{script_dir}/submit_gau.sh", "-i", f"{basename}.com"]
    log_opted = False
    if os.path.isfile(f"{basename}.log"):
        if "Normal termination" in open(f"{basename}.log", 'r').readlines()[-1]:
            log_opted = True
            logging.info(f"gau opt for {basename} is done")
            results = subprocess.run(["echo","opt done"],capture_output=True,text=True)
    if not log_opted:
        jobid = subprocess.run(_cmd, capture_output=True, text=True).stdout.strip().split()[-1]
        time.sleep(5)
        #query = f"sacct -X --parsable2 --noheader --format JobID,State%30,Elapsed,End -j {jobid}"
        query = f"squeue -j {jobid}"
        while True:
            time.sleep(30)
            results = subprocess.run(query.split(), capture_output=True, text=True)
            #if "RUNNING" not in results.stdout: break
            if results.returncode != 0: break # means that job is not in active queue!
    logging.debug(f"gaussian results.stderr:\n{results.stderr}")
    return f"{basename}.log"

def handle_linear_angle(logFile,SD):
    _cmd = [f"{SD}/extract_last_pdb.sh", logFile, "LIG"]
    results = subprocess.run(_cmd, capture_output=True, text=True)
    pdbFile = logFile[:-4]+"_linear_angle.pdb"
    with open(pdbFile, "w") as fp:
        fp.write(results.stdout)
    return pdbFile

def calc_esp(args):
    # 3 opts: opt=None, opt=gau, opt=psi4
    # 2 esps: esp=gau, esp=psi4
    pdb4psi = args.inp_file
    basename = pdb4psi[:-4]
    espFile = f"results_{basename}.esp"
    pdbFile = pdb4psi
    if not os.path.isfile(espFile):
        if args.optimization in ["None","none"]:
            pdb4psi = args.inp_file
            if args.esp_engine == "psi4":
                espFile = psi4_esp_run(inp_file=pdb4psi,cpus=args.cpus,
                lig_charge=args.lig_charge,optimization=0,mems=args.mems,
                functional=args.esp_functional,basis=args.esp_basis)
            elif args.esp_engine == "gau":
                gen_gau_inp(pdb_file=pdb4psi,charge=args.lig_charge,multip="1",
                    func=args.esp_functional,basis=args.esp_basis,
                    pre_opt=0,just_opt=0,cpus=args.cpus,mems=args.mems)
                logFile = run_gau(basename,args.cpus)
                logging.debug(f"gaussian finished: logFile is {logFile}")
                pdbFile = pdb4psi
                espFile = f"results_{basename}.esp"
        elif args.optimization == "gau":
            if args.esp_engine == "psi4":
                gauOpts = ["B3LYP","6-31G**",0,1]
            elif args.esp_engine == "gau":
                gauOpts = [args.esp_functional,args.esp_basis,1,0]
            gen_gau_inp(pdb_file=pdb4psi,charge=args.lig_charge,multip="1",
                func=gauOpts[0],basis=gauOpts[1],cpus=args.cpus,mems=args.mems,
                pre_opt=gauOpts[2],just_opt=gauOpts[3])
            logFile = run_gau(basename,args.cpus)
            logging.debug(f"gaussian finished: logFile is {logFile}; pdb4psi[:-4]: {pdb4psi[:-4]}")
            print(f"Check: path: {os.getcwd()}, logFile: {logFile}, espFile: {espFile}",flush=True)
            if "Bend failed for angle" in open(logFile, 'r').read():
                logging.debug(f"{logFile} is not properly optimized due to linear angle, re-sub one more time!")
                tmpPdbF = handle_linear_angle(logFile,args.script_dir)
                gen_gau_inp(pdb_file=tmpPdbF,charge=args.lig_charge,multip="1",
                            func=gauOpts[0], basis=gauOpts[1], cpus=args.cpus, mems=args.mems,
                            pre_opt=gauOpts[2], just_opt=gauOpts[3]
                            )
                logFile = run_gau(tmpPdbF[0:-4],args.cpus)

            #_cmd = [f"{args.script_dir}/extract_last_pdb.sh", f"{pdb4psi[:-4]}.log", f"{args.resName}"]
            _cmd = [f"{args.script_dir}/extract_last_pdb.sh", logFile, f"{args.resName}"]
            results = subprocess.run(_cmd,capture_output=True,text=True)
            logging.info(f"gaussian opt done: stderr: {results.stderr}")
            os.system(f"{args.script_dir}/extract_last_xyz.sh  {logFile} > {basename}.xyz")
            pdbFile = f"{basename}_gauopt.pdb"
            with open(pdbFile,"w") as fp:
                fp.write(results.stdout)
            logging.info(f"stderr of extracting last geom from gaussian log: {results.stderr}")
            opted_pdb = f"{basename}_opt"
            with open(f"{opted_pdb}.pdb", 'w') as fp:
                fp.write(''.join(mapping_atom_back(f"{basename}.xyz", f"{basename}.pdb")))
            if args.just_opt: return [f"{opted_pdb}.pdb"]
            if args.esp_engine == "gau":
                espFile = f"results_{basename}.esp"
            elif args.esp_engine == "psi4":
                espFile = psi4_esp_run(inp_file=pdbFile,cpus=args.cpus,
                lig_charge=args.lig_charge,optimization=0,mems=args.mems,
                functional=args.esp_functional,basis=args.esp_basis)
        elif args.optimization == "psi4":
            pdbFile = f"{basename}_opt.pdb"
            espFile = psi4_esp_run(inp_file=pdbFile,
                lig_charge=args.lig_charge,optimization=1,cpus=args.cpus,
                functional=args.esp_functional,basis=args.esp_basis,mems=args.mems)

    logging.info(f"psi4_esp_run results files: {espFile}")
    if os.path.isfile(espFile):
        if not os.path.isfile(f"{espFile[:-4]}.esp.pdb"):
            gesp2pdb = subprocess.run([f"{args.script_dir}/gesp2pdb.py",pdbFile,espFile],capture_output=True)
            print(f"gesp2pdb results: {gesp2pdb.stdout}\n{gesp2pdb.stderr}")
        cmd.load(f"{espFile[:-4]}.lig.pdb")
        cmd.load(f"{espFile[:-4]}.esp.pdb")
        cmd.set("sphere_transparency", 0.5)
        cmd.set("sphere_scale", 0.1)
        cmd.spectrum("b", "red_white_blue", f"{espFile[:-4]}.esp")
        cmd.label(f"{espFile[:-4]}.lig", "'%1.2f'%b")
        cmd.set("label_font_id", 1)
        cmd.set("label_size", 20)
        cmd.set("bg_rgb","white")
        cmd.set("stick_h_scale", 1)
        cmd.show("sticks",f"{espFile[:-4]}.lig")
        cmd.save("esp_charge_distribution.pse")
        cmd.reinitialize()
    return [pdbFile,espFile]

def process_covalent_pose(args):
    # ["covComPdb4tleap.pdb","capped_cys.pdb",
    # "lig_retained.pdb","automated_lig_capped.pdb"]
    pdbList = prepare_capp(pdb=args.inp_file,ligSele=args.ligSel)
    pdb4psi = pdbList[-1]
    for line in open(pdbList[-2],'r'): # lig_retained.pdb
        if line[0:6] in ['ATOM  ','HETATM']:
            args.resName = line[17:20].strip()
            break
    args.inp_file = pdb4psi
    logging.debug(f"Checking files: {pdbList}")
    if args.charge_type == "bcc":
        logging.warning("BCC charge is not well supported for covalent ligand")
    if args.inp_file[-3:] != "pdb":
        logging.error("only complex file in PDB format is supported for -pt cov")
    if args.optimization == "psi4" and args.esp_engine == "gau":
        logging.error("optimization by psi4 and esp by gaussian is not supported")
    elif args.optimization == "psi4":
        logging.warning("Optimization by psi4 is slow and difficult in convergence")

    if args.charge_type.lower() == "esp":
        # fitting in one command
        files = calc_esp(args)
        return files
    elif args.charge_type.lower() == "resp":
        files = calc_esp(args)
        _cmd = [f"{args.script_dir}/fit_cov_resp.py",pdb4psi,files[1],args.resName]
        logging.info(f"running resp fitting: {_cmd}")
        results = subprocess.run(_cmd,capture_output=True,text=True)
        logging.info(f"resp fitting stdout: {results.stdout}")
        return results
    elif args.charge_type.lower() == "bcc":
        #LL's flavor, LIS as a fake protein residue
        _cmd = [f"{args.script_dir}/fit_cov_oneReside.py",pdb4psi,args.lig_charge,"bcc"]
        logging.info(f"running: {' '.join(_cmd)} now, please do not perform opt!")
        # results = ["cappedPdb.pdb", oldNames, newNames]
        results = subprocess.run(_cmd,capture_output=True,text=True)
        #TODO: if optimization performed, recover the atom names back in prepi!
        return results

def process_normal_pose(args):
    if args.ligSel == "organic":
        logging.warning("ligSel is not given, could be problematic in default settings")
    cmd.load(args.inp_file)
    cmd.select("ligSel", args.ligSel)
    resSet = set()
    atomList = cmd.get_model("ligSel").atom
    res_List = [resSet.add((x.chain,x.resn,x.resi)) for x in atomList]
    if len(resSet) > 1:
        logging.error(f"more than 1 ligands fould in the {args.ligSel}"+
                      f"please explicity specify the ligand selection: e.g. A:LIG:442")

    pdbFile = f"{args.inp_file[:-4]}_{atomList[0].resn}.pdb"
    cmd.set("retain_order",1)
    cmd.save(pdbFile,args.ligSel)
    cmd.reinitialize()
    args.inp_file = pdbFile
    process_single_ligand(args)

def process_single_ligand(args):
    basename = args.inp_file[:-4]
    if args.charge_type == "bcc":
        if args.inp_file[-3:] != "mol2":
            cmd.load(args.inp_file)
            args.resName = cmd.get_model("all").atom[0].resn
            cmd.set("retain_order", 1)
            cmd.save(f"{basename}.mol2")
            cmd.reinitialize()
        acFlags = args.ac_flags.replace("@","-")
        if "dr no" in acFlags:
            cmds = [
                f'''antechamber -i {basename}.mol2 -fi mol2 -o {basename}.ac -fo ac  \
                            -c bcc -s 2 -at gaff -nc {args.lig_charge} -rn {args.resName}  -pf Y {acFlags}''']
        else:
            cmds = [
                f'''antechamber -i {basename}.mol2 -fi mol2 -o {basename}.ac -fo ac  \
                            -c bcc -s 2 -at gaff -nc {args.lig_charge} -rn {args.resName}  -pf Y {acFlags}''',
            f"prepgen -i {basename}.ac -o {basename}.prepi -f inp -rn {args.resName}",
            f"parmchk2 -i {basename}.prepi -o {basename}.frcmod -f prepc -a Y"
        ]
        for _cmd in cmds:
            logging.info(f"running cmds for bcc charge fitting:")
            logging.info(_cmd)
            results = subprocess.run(_cmd.split(),capture_output=True,text=True)
            logging.debug(f"results: {results.stdout}\n{results.stderr}")
        if os.path.isfile("sqm.pdb"):
            if len(open("sqm.pdb","r").readlines()) > 2: # should have at least 2 atoms
                cmd.load('sqm.pdb')
                cmd.set("retain_order", 1)
                cmd.save('sqm.xyz')
                cmd.reinitialize()
        files = ['sqm.xyz']
    else:
        if args.inp_file[-3:] != "pdb":
            cmd.load(args.inp_file)
            args.resName = cmd.get_model("all").atom[0].resn
            cmd.set("retain_order", 1)
            cmd.save(f"{basename}.pdb")
            cmd.reinitialize()
        args.inp_file = f"{basename}.pdb"
        files = calc_esp(args)
        if args.just_opt: return files
        if args.charge_type == "resp":
            if args.inp_file[-3:] == "pdb": # recommand PDB format for RESP!
                files[0] = args.inp_file
            pseOption = " "
            if os.path.isfile("esp_charge_distribution.pse"):
                pseOption = "-pse esp_charge_distribution.pse"
            _cmd = f'''{args.script_dir}/fitresp_from_esp.sh -pre {basename} 
                   -rn {args.resName} -e {files[1]} -i {files[0]} -pf {pseOption}'''
            results = subprocess.run(_cmd.split(),capture_output=True,text=True)
            logging.info(f"_cmd: {_cmd}\nstdout: {results.stdout}\nstderr: {results.stderr}")
            return results
        elif args.charge_type == "esp":
            # maybe bugs here, need to be careful!
            if os.path.isfile("esp_charge_distribution.pse"):
                return files
    return files

def main():
# Main starts here:
    args = read_par()
    args = args.parse_args()
    logging.info(f"working in {args.workDir}")
    logging.info(f"script_dir is {args.script_dir}")
    os.makedirs(args.workDir,exist_ok=True)
    if os.path.isfile(args.inp_file) or os.path.islink(args.inp_file):
        cp_cmd = ["/usr/bin/cp","-vf",args.inp_file,args.workDir]
        res = subprocess.run(cp_cmd,capture_output=False,shell=False)
        logging.info(f"copying inp to workdir:\n{res.stdout}\n{res.stderr}")
    os.chdir(args.workDir)
    if os.path.isfile(args.inp_file) or os.path.islink(args.inp_file):
        if args.pdb_type == "cov":
            process_covalent_pose(args)
        elif args.pdb_type == "com":
            process_normal_pose(args)
        elif args.pdb_type == "lig":
            process_single_ligand(args)
    else:
        sdfFile = smi2sdf(args.inp_file)
        args.inp_file = sdfFile
        process_single_ligand(args)

if __name__ == "__main__":
    main()
