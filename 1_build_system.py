#!/public/home/lijunhao/soft/conda_envs/fldev/bin/python
import glob
import os
import subprocess
import sys
import json
import time
import socket
from pathlib import Path

script_dir = Path(__file__).resolve().parent
print(script_dir)

# need to adjust the module path
if  '192.168' in socket.gethostname():
    module_py = '/opt/homebrew/Cellar/lmod/9.2/init/env_modules_python.py'
    ambertool = 'ambertools/24.0'
else:
    module_py = '/public/software/apps/module/init/python.py'
    ambertool = 'amber/a22t23'
exec(open(module_py).read())
module('purge') # type: ignore
module(f'load {ambertool}') #type: ignore
module('list') # type: ignore
# if you don't have Module/lmod, comment the about codes
# and make sure Ambertools is globally installed

def tleap_PP(receptor="full_st_nh_leap_ter.pdb",ff='ff14SB',leap_file=[],
             salt_num=0,wat_box='TIP3P',ligFF='gaff'):
    pre_loading = pre_loads(leap_file)
    bondInfo = "#pymol exported structure with CONECT record, no bond text"
    if salt_num == 0:
        salt_text = "  "
    else:
        salt_text = f'''
        addionsrand cpx Na+ {salt_num}
        addionsrand cpx Cl- {salt_num}
        '''
    txt = f'''
    # non-standard DNA/RNA residue support is ongoing...
    source leaprc.protein.{ff}
    source leaprc.RNA.OL3
    source leaprc.DNA.OL15
    source leaprc.phosaa14SB
    source leaprc.water.{wat_box.lower()}
    source leaprc.{ligFF}
    {pre_loading}
    cpx = loadpdb {receptor}
    {bondInfo}
    savepdb cpx cpx_vac.pdb
    solvatebox  cpx {wat_box.upper()}BOX 12
    addionsrand cpx Na+ 0
    addionsrand cpx Cl- 0
    {salt_text}
    saveamberparm cpx cpx_solvated.prmtop cpx_solvated.inpcrd
    savepdb cpx cpx_solvated.pdb
    quit
    '''
    check = f'''
    source leaprc.protein.{ff}
    source leaprc.RNA.OL3
    #source leaprc.DNA.bsc1
    source leaprc.DNA.OL15
    source leaprc.phosaa14SB
    source leaprc.water.{wat_box.lower()}
    source leaprc.{ligFF}
    {pre_loading}
    cpx = loadpdb {receptor}
    {bondInfo}
    savepdb cpx cpx_vac.pdb
    check cpx
    quit
    '''
    with open('tleap_PP.in', 'w') as f:
        f.write(txt)
    with open('check_cpx.in', 'w') as fp:
        fp.write(check)
    return ['check_cpx.in','tleap_PP.in']

def pre_loads(leap_file):
    load_text =[]
    for f in leap_file:
        if '.prep' in f:
            load_text.append(f'loadamberprep {f}\n')
        if '.lib' in f:
            load_text.append(f'loadoff {f}\n')
        if '.frcmod' in f:
            load_text.append(f'loadamberparams {f}\n')
    return ''.join(load_text)

def define_esp_mem(pdb):
    n_atom = 0
    memory = 80
    for line in open(pdb,'r'):
        if line[0:6] in ['ATOM  ','HETATM']:
            n_atom += 1
    if n_atom < 40: memory = 80
    elif n_atom < 80: memory = 120
    elif n_atom < 120: memory = 160
    return [n_atom,memory]

def get_num_ions(log,conc):
    for line in open(log,'r'):
        if 'Volume: ' in line:
            y = float(line.split()[1])
            x = float(conc)*1000
            return int((y/(10**27))*(x/1000*6.022*(10**23)))

def main(pdb_name,lig_info,lig_charge='resp',ff="ff14SB",
    wat_box="TIP3P",salt_conc="0"):
    # lig_info is deprecated now...
    print(f"capping and charge guessing for blocks: {pdb_name}:",flush=True)
    dict_cmd = f"{script_dir}/2_prepare_pep_caps.py {pdb_name}"
    print(dict_cmd,flush=True)
    resp_dict_run = subprocess.run(dict_cmd.split(),capture_output=True,text=True)
    print(f"stdout:\n{resp_dict_run.stdout}\nstderr:\n{resp_dict_run.stderr}",flush=True)
    resp_dict = json.load(open("blocks_charges.json", 'r', encoding="utf-8"))
    print(f"resp_dict is: {resp_dict}")
    leap_file, job_ids = [], []
    for f in resp_dict.keys():
        #chg = lig_info[f] # passed
        chg = str(resp_dict[f])
        pre = f[0:-4]
        if f.startswith('pure_lig'):
            #ct = 'resp'# time to add support for bcc
            ct = lig_charge.lower()
            resn = pre.split('_')[-1]
        else: # NCAAs
            ct = 'esp' if lig_charge == 'resp' else lig_charge
            resn = 'LIG' # only fits esp, doesn't matter what resn it is
        _mem = str(define_esp_mem(f)[1])
        wd = f"{lig_charge}_{pre}"
        # one-by-one is too slow!
        #cmd = [f"{scriptDir}/toolkits/tool2_charge_fitting/atomic_charge_fitting.py",
        cmd = [f"{script_dir}/3_PPmd_gp4.sh", os.getcwd(),
               "-wd",wd,"-pt","lig","-inp",f,"-opt","gau","-chg",chg,
               "-ct",ct,"--resn",resn,"--cpus","24","--mems",_mem]
        espout = glob.glob(f"{wd}/result*.esp")
        if len(espout) > 0 and os.path.isfile(espout[0]):
            if 'ESP VALUES AND GRID' in open(espout[0]).read():
                print(f"ESP file found: {espout[0]}\nIf error occurs, remove the old esp to re-run!",flush=True)
                continue
        print(f"### running cmd: {' '.join(cmd)}")
        res = subprocess.run(cmd,capture_output=True,text=True)
        job_ids.append(res.stdout.strip().split()[-1])
    time.sleep(5)
    print(f"esp/resp fitting jobs submitted, jobIds: {job_ids}")
    while len(job_ids) > 0:
        for jid in job_ids:
            #sacct_run = f"sacct -X --parsable2 --noheader --format State -j {jid}"
            squ_run = f"squeue -j {jid}"
            squ_out = subprocess.run(squ_run.split(),capture_output=True,text=True)
            # could be failed, check carefully of the work dir when the job failed!!!
            #if sacct_out.stdout.strip() not in ['RUNNING','PENDING']:
            if squ_out.returncode != 0: # not in the active job queue!
                job_ids.remove(jid)
        time.sleep(10)
    print(f"looks like that all the charge fitting slurm jobs ended")
    time.sleep(10)

    for f in resp_dict.keys(): #standard capping
        pre = f[0:-4]
        wd = f"{lig_charge.lower()}_{pre}"
        if f.startswith('pure_lig'):
            leap_file += glob.glob(f"{wd}*/{pre}*prep*")
            leap_file += glob.glob(f"{wd}*/{pre}*frcmod")
            continue
        os.makedirs(wd,exist_ok=True)
        os.chdir(wd)
        print(f"### preparing parameters for {f} in {os.getcwd()}")
        if os.path.isfile(f"{pre}.log"):
            os.system(f"{script_dir}/extract_last_xyz.sh {pre}.log > {pre}.xyz")
        if os.path.isfile(f"{pre}_opt_init.xyz"):
            os.system(f"cp -v {pre}_opt_init.xyz {pre}.xyz")
        #_cmd = [f"{scriptDir}/md/P_P_MD/fit_pep_chgs.py",f,f"results_{pre}*.esp", f"{pre}.xyz"]
        if lig_charge == 'resp':
            espF = glob.glob(f"results_{pre}*.esp")
            _cmd = f"{script_dir}/4_fit_pep_chgs.py -pdb {f} -esp {espF[0]} -xyz {pre}.xyz"
        elif lig_charge == 'bcc':
            ### BCC charge is still not working well, use it carefully!
            _cmd = f"{script_dir}/4_fit_pep_chgs.py -pdb {f} -chg_type bcc -xyz {pre}.xyz"
        print(f"### running cmd: {_cmd}")
        res = subprocess.run(_cmd.split(),capture_output=True,text=True)
        leap_file += [f"{wd}/{x}" for x in res.stdout.strip().split("\n")[-1].split()]
        print(f"response: {res.stdout}\n{res.stderr}")
        os.chdir("..")

    lig_ff = 'gaff'
    #prep = pdb_pre_process(pdb_name)
    #generate_lig_top(lig_ff,lig_charge,lig_info)
    recp_pdb = "full_st_nh_leap_ter.pdb" # output from 2_prepare_pep_caps.py
    pdb4amb = f"pdb4amber -i {recp_pdb} -o p4a.pdb -l p4a.log" # maybe CYX needs to present at the beggining!
    res = subprocess.run(pdb4amb.split(),capture_output=True,text=True)
    #recp_pdb = 'p4a.pdb'
    print(f"response of pdb4amb: {res.stdout}\n {res.stderr}",flush=True)
    print(f"files list loading before cpx: {leap_file}",flush=True)
    leapF = tleap_PP(recp_pdb,ff,leap_file,salt_num=0,wat_box=wat_box,ligFF=lig_ff)
    cyclic = f"{script_dir}/cyclic_tleap.sh {leapF[0]} {recp_pdb} {leapF[1]}"
    print(f"running cmd: {cyclic}",flush=True)
    os.system(cyclic)
    subprocess.run(f"tleap -s -f {leapF[1]}".split(),capture_output=True,text=True)
    if float(salt_conc) > 0:
        ion_num = get_num_ions('leap.log',salt_conc)
        tleap_PP(recp_pdb,ff,leap_file,salt_num=ion_num,wat_box=wat_box,ligFF=lig_ff) # type: ignore
        os.system(cyclic)
        subprocess.run(f"tleap -s -f {leapF[1]}".split(),capture_output=True,text=True)

if __name__ == "__main__":
    #lig_info e.g. "ZDC:0,DCY:0,DTR:0,DLY:1,8VH:0"
    pdb_name: str = sys.argv[1]
    lig_info=sys.argv[2].split(',')
    lig_charge=sys.argv[3].lower() # resp
    ff=sys.argv[4]
    temp = sys.argv[5]
    water_box = sys.argv[6]
    job_time=sys.argv[7] # unit is hour
    salt_conc = sys.argv[8]
    #md_in_set = sys.argv[9]
    main(pdb_name,lig_info,lig_charge,ff,water_box,salt_conc)
