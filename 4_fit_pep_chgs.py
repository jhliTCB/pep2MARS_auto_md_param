#!/public/home/lijunhao/soft/conda_envs/fldev/bin/python
# The pdb file feeded to this step should have unified residue names and atom names!
# The input must be a prepared complex with all hydrogen atoms!
# for the pari, triple, and quartet, need to know if there are ligand inside!
import argparse
import json
import os
import subprocess
import math
import sys
#sys.path.append('/public/home/lijunhao/soft/pymol_schrd3151/lib/python3.10/site-packages')
from pymol import cmd, stored
#from prepare_pep_caps import guessing_charge

# for debug
exec(open('/public/software/apps/module/init/python.py').read())
module('purge')
module('load amber/a22t23 gaussian/16.c02')
module('list')

global ACE, NME, AMBTYPE
# covering more atom names from different source
ACE = {'CH3':-0.3662,'H1':0.1123,'H2':0.1123,'H3':0.1123,
       'C':0.5972,'O':-0.5679,'CA':-0.3662,'HA':0.1123,
       'HA1':0.1123,'HA2':0.1123, 'HA3':0.1123, '1HH3':0.1123,
       'H01':0.1123,'H02':0.1123, 'H03':0.1123, 'H04':0.1123,
       '2HH3':0.1123, '3HH3':0.1123}
NME = {'N':-0.4157,'H':0.2719,'HNN':0.2719,'C':-0.1490,'CA':-0.1490,
       'H1':0.0976,'H2':0.0976,'H3':0.0976,'HA':0.0976,
       'H01':0.0976,'H02':0.0976,'H03':0.0976,'H04':0.0976,
       'HA1':0.0976,'HA2':0.0976,'HA3':0.0976, 'CH3':-0.1490,
       '1HH3':0.0976, '2HH3':0.0976, '3HH3':0.0976}
# only mapping the BB atomtypes, for other DUs, make gaff type upper!
AMBTYPE = {'N':'N','H':'H','C':'C','O':'O','CA':'CX'}

def repalce_mol2_DU(inp,full_ac,indices):
    # single residue itself must be open-valence, thus DU occurs!
    # the backbone of NCAAs sometimes are not correctly recognized
    # inp and gaff have same number of lines: two mol2 from same pdb
    ante = "%7d %-8s %10.4lf %10.4lf %10.4lf %-6s %5d %-8s %9.6lf\n"
    flag = False
    ac_info = get_ac_info(full_ac)
    out,ind = [],0
    for line in open(inp,'r'):
        if 'TRIPOS>ATOM' in line:   flag = True
        elif 'TRIPOS>BOND' in line: flag = False
        if flag and 'TRIPOS>ATOM' not in line:
            atm = line.split()
            # prone to out-range if atom order changed, make sure antechamber not change atom order for ac file!
            new_type = ac_info[indices[ind]][-1]
            line = ante % (int(atm[0]),atm[1],float(atm[2]),float(atm[3]),
                float(atm[4]), new_type, int(atm[6]), atm[7],float(atm[8]))
            ind += 1
        out.append(line)
    return out

def get_init_chg(pdb):
    # return a list masking the positions of capped atoms
    # read the REMARK of residue-residue covalent bonds added by pymol!
    # and split the pdb into separated units!
    # ligResName is removed, because ligand/linker should not have NHCO names
    NHCO = ['N', 'H', 'C', 'O']
    ind, n_cap = 0, 0
    resp0, bondinfo, keyinfo = [], [], []
    ressInfo, ressType = {}, {} # true non-capped atoms
    #resn_list = get_resn_list(pdb)
    with open(pdb, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if 'REMARK: bond' in line: bondinfo.append(line.split(':')[-1])
        if line[0:6] in ['ATOM  ', 'HETATM']:
            atmNam = line[12:16].strip()
            resNam = line[17:20].strip()
            keyinfo.append([atmNam,resNam,ind+1])
            #resNum = line[22:26].strip()
            #chainN = line[21:22]
            if resNam == 'ACE':
                resp0.append(ACE['H1'] if atmNam[0] in ['H', 'h'] else ACE[atmNam])
                n_cap += 1
                #cap.append(replace_cap_mc(line,atmNam))
            elif resNam in ['NME','NMA']:
                resp0.append(NME[atmNam])
                #cap.append(replace_cap_mc(line,atmNam))
                n_cap += 1
            else:
                resp0.append(0.0000) # resp0 for resp1 fitting
                if resNam not in ressInfo.keys():
                    ressInfo[resNam] = []
                ressInfo[resNam].append([ind,line])
            ind += 1

    for res in ressInfo.keys():
        nhco = [x for x in ressInfo[res] if x[1][12:16].strip() in NHCO]
        if len(nhco) == 0:
            ressType[res] = 'gaff'
        else:
            ressType[res] = 'amber'
    print("number of capped atoms are {}".format(n_cap))
    return [resp0,ressInfo,ressType,bondinfo,keyinfo]

def write_init_chg(chgi, fout='qin.tmpl'):
    fp = open(fout, 'w')
    for i in range(len(chgi)):
        j = i + 1
        fp.write("{:10.6f}".format(float(chgi[i])))
        if float("{:.1f}".format(math.fmod(j,8))) == 0.0:
            fp.write('\n')
    fp.close()
    return 0

def mod_resp_inp(chgi,respin,stage=1,fout='resp1.in'):
    # Add a '-99' constraint to the capped atoms, because their
    # charges are fixed and summed to zero
    i, j = -2, 0
    flag = False
    fp = open(fout,'w')
    for line in open(respin,'r'):
        j += 1
        if "Resp charges for" in line and j > 1:
            flag = True
            i += 1
            fp.write(line)
            continue
        if not flag:
            fp.write(line)
            if 'qwt' in line and stage == 1:
                fp.write(' iqopt = 2,\n')
        elif i == -2:
            fp.write(line)
        elif i == -1:
            fp.write(line)
            i += 1
            continue
        elif len(line.split()) > 1:
            #print(i)
            atmic = int(line.split()[0])
            state = line.split()[1]
            if chgi[i] == 0.0:
                fp.write(line)
            else:
                fp.write('{:5d}{:5d}\n'.format(atmic,-99))
            i += 1
        else:
            fp.write(line)
    fp.close()
    return 0

def replace_amb_du_type(gaff_ac, amber_ac):
    '''even from a correct .esp file, using amber atom type by antechamber can still generate DU atom type,
    just repalce them with the corresponding one in upper case from gaff type'''
    new_ac = []
    for line1, line2 in zip(open(gaff_ac,'r'),open(amber_ac,'r')):
        if 'DU' in line2.split()[-1]:
            line2 = line2.replace('DU\n',f'{line1.split()[-1].upper()}\n')
        new_ac.append(line2)
    with open(amber_ac,'w') as fp: fp.write(''.join(new_ac))

def get_init_all_resp(esp,init_chg):
    base_esp = '.'.join(esp.split('.')[0:-1])
    write_init_chg(init_chg, fout='cap.chg')
    print(f"init_chg info for 2-step resp: {init_chg}",flush=True)
    print(f"sum of caps charges={sum(init_chg)},len(init_chg)={len(init_chg)}",flush=True)

    os.makedirs('tmp', exist_ok=True)
    cmdx = f"$AMBERHOME/bin/espgen -i {esp} -o {base_esp}.dat"
    print(cmdx,flush=True)
    os.system(cmdx)
    os.chdir('tmp')
    print(f"now working in {os.getcwd()}",flush=True)
    for atype in ['gaff','amber']:
        cmd0 = f"antechamber -seq no -an Y -fi gesp -i ../{esp} -fo ac -o full_capped_resp_{atype}.ac -c resp -at {atype}"
        print(cmd0,flush=True)
        os.system(cmd0)
    replace_amb_du_type('full_capped_resp_gaff.ac','full_capped_resp_amber.ac')
    if os.path.isfile('ANTECHAMBER.ESP'):
        eDat = 'tmp/ANTECHAMBER.ESP'
    else:
        eDat = f"{base_esp}.dat"
    os.chdir('..')
    print(f"now working in {os.getcwd()}", flush=True)
    mod_resp_inp(init_chg, 'tmp/ANTECHAMBER_RESP1.IN', stage=1, fout='resp1.in')
    mod_resp_inp(init_chg, 'tmp/ANTECHAMBER_RESP2.IN', stage=2, fout='resp2.in')
    cmd1 = f"$AMBERHOME/bin/resp -O -i resp1.in -o resp1.out -p resp1.pch -t resp1.chg -q cap.chg -e {eDat}"
    cmd2 = f"$AMBERHOME/bin/resp -O -i resp2.in -o resp2.out -p resp2.pch -t resp2.chg -q resp1.chg -e {eDat}"
    print(f"running customized resp fit\n{cmd1}\n{cmd2}\n",flush=True)
    os.system(cmd1)
    os.system(cmd2)
    # now the resp2.chg can be used, the sum charges of the capped atoms should zero.
    all_resp = []
    for line in open('resp2.chg', 'r'):
        all_resp += [float(x) for x in line.strip().split()]
    return all_resp

def get_ac_info(acF):
    ac_info = []
    for line in open(acF, 'r'):
        if 'ATOM' in line:
            # ac.c: fprintf(fpout, "ATOM%7d  %-4s%-4s%5d%12.3f%8.3f%8.3f%10.6lf%10s\n"
            atom_name = line[13:17].strip()
            resname = line[17:21].strip()
            atom_id = line.split()[1] # safe, won't have atoms more than 10^6 !
            atom_type = line.split()[-1] # safe because no type name can longer than 10 chars!
            ac_info.append([atom_name,resname,atom_id,atom_type])
    return ac_info

def label_cap_for_bcc(pdb):
    '''A better way is to label tmp_full_block.ac
    but it will need to handle the atom names in connection records!'''
    out = []
    for line in open(pdb, 'r'):
        if line[0:6] in ['ATOM  ', 'HETATM']:
            if 'NME' in line or 'ACE' in line or 'NMA' in line:
                old_name = line[12:16]
                old_repl = line[12:26]
                if len(old_name.strip()) == 4:
                    new_name = f"_{old_name[1:]}"
                else:
                    new_name = "{:^4s}".format("_"+old_name.strip())
                new_repl = old_repl.replace(old_name,new_name)
                line = line.replace(old_repl,new_repl)
        out.append(line)
    return out

def get_prepi_info(prepiF):
    ''' prepgen.c: # AmberTools >= 22
     if (cartindex == 0) {
            fprintf(fprep, "%4d  %-5s %-5s %c", j + 1, newatom[j].name,
                    newatom[j].ambername, aatype[j].name[0]);
            fprintf(fprep, "%5d%4d%4d%10.3f%10.3f%10.3f%10.6f\n", bond + 1, angle + 1,
                    twist + 1, bondv, anglev, twistv, newatom[j].charge);
        }
    '''
    prepi_info = []
    for line in open(prepiF, 'r'):
        if len(line) > 70 and 'DUMM' not in line:
            atom_name = line.split()[1]
            atom_type = line.split()[2]
            partial_charge = line[62:72].strip()
            prepi_info.append([atom_name,partial_charge,atom_type])
    return prepi_info

def get_init_all_bcc(pdb,sdf,resn,net_chgs,key_info):
    '''TODO: It is not working for blocks with more than two caps! need to replace it with cpptraj'''
    basename = '.'.join(pdb.split('.')[0:-1])
    stored.resn = []
    cmd.load(sdf)
    cmd.set("retain_order",1)
    cmd.save(f'{basename}.mol2')

    cmd.iterate("enabled","stored.resn.append(resn)")
    ress = sorted(set(stored.resn))
    print(f"resnames in ress: {ress}\n",flush=True)
    res_chg_list = []
    for res in ress:
        if res in ['ACE','NME']: continue
        res_atm_lst = [x.name for x in cmd.get_model(f"resn {res}").atom]
        cmd.create(f"{res}_e1",f"r. {res} extend 1")
        res_ext_lst = [x for x in cmd.get_model(f"{res}_e1 & ! r. {res}")]
        omit_names = []
        for atm in res_ext_lst:
            cmd.edit(f"r. {atm.resn} & r. {atm.name}")
            cmd.h_fill()
            cmd.edit()
            for ext in [x for x in cmd.get_model(f"r. {atm.resn} & r. {atm.name}")]:
                if ext.name in res_atm_lst or ext.name in omit_names:
                    for i in range(100):
                        newN = ext.name[0]+str(i)
                        if newN not in res_atm_lst and newN not in omit_names:
                            omit_names.append(newN)
                            cmd.alter(f"r. {ext.resn} & r. {ext.name}",f"name='{newN}'")
                            break
            cmd.alter(f"r. {atm.resn} & {res}_e1",f"resn='{res}")
        cmd.set("retain_order",1)
        cmd.save(f"{res}_e1.mol2",f"{res}_e1")
        cmd.save(f"{res}_e1.pdb",f"{res}_e1")
        res_chg_list.append(f'{res}_e1.pdb')
        with open(f"{res}_e1.mc", 'w') as fp:
            fp.write('\n'.join([f"OMIT_NAME {x}" for x in omit_names]))
            fp.write('CHARGE to_be_seded\n')

    with open(f'_{basename}.pdb','w') as fp:
        fp.write(''.join(label_cap_for_bcc(f'{basename}.pdb')))
    cmd.load(f'_{basename}.pdb') # pdb file without any residue information
    cmd.set("retain_order", 1)
    cmd.save(f'_{basename}.mol2')
    cmd.reinitialize()
    tmp_ac = "tmp_full_block.ac"
    _cmd = f"antechamber -fi mol2 -i _{basename}.mol2 -fo ac -o {tmp_ac} -nc {net_chgs} -c bcc -dr no"
    print(_cmd,flush=True)
    if not os.path.isfile(tmp_ac): os.system(_cmd)
    if not os.path.isfile(tmp_ac):
        print(f"error in generating a full ac file!",flush=True)
        exit(1)
    ac_info = get_ac_info(tmp_ac)
    cap_ind = [i for i in range(len(key_info)) if key_info[i][1] in ['ACE','NME']]
    #omit_name = []
    #for oo in cap_ind:
    #    if ac_info[oo][0] not in omit_name:
    #        omit_name.append(ac_info[oo][0])
    omit_name = [ac_info[i][0] for i in cap_ind]
    with open('omit_name.mc','w') as omit_file:
        omit_file.write('\n'.join([f"OMIT_NAME {x}" for x in omit_name]))
        omit_file.write(f'\nCHARGE {net_chgs}\n')
    _cmd = "prepgen -i tmp_full_block.ac -o _decapped_.prepi -m omit_name.mc -rn TMP"
    print(f"going to run cmd {_cmd} to get _decapped_prepi",flush=True)
    os.system(_cmd)
    prepi_info = get_prepi_info('_decapped_.prepi')
    print(f"prepi info is {prepi_info}",flush=True)
    all_bcc = []
    for i in range(len(key_info)):
        if key_info[i][1] in ['ACE','NME']: all_bcc.append('0.000000')
        tmp_name = ac_info[i][0]
        for j in range(len(prepi_info)):
            if prepi_info[j][0] == tmp_name:
                all_bcc.append(prepi_info[j][1])
    print(f"\n_decapped_.prepi charges order: {[x[1] for x in prepi_info]}\n",flush=True)
    print(f"sum of charges in _decapped_.prepi: {sum([float(x[1]) for x in prepi_info])}\n",flush=True)
    print(f"ordered charges from original PDB: {all_bcc}\n",flush=True)
    print(f"sum of charges after re-order: {sum([float(x) for x in all_bcc])}\n",flush=True)
    return all_bcc

def gen_right_ac(ac1,ac2,indicies):
    '''ac1: the incorrect lig.ac, ac2: the right gaff ac file in tmp dir'''
    ac2_info = get_ac_info(ac2)
    ind, out = 0, []
    for line in open(ac1,'r'):
        if line.startswith('ATOM'):
            type_str = line[70:]
            full_ind = indicies[ind]
            new_type = f"{str(ac2_info[full_ind][-1]):>{len(type_str)}}"
            line = line[0:70] + new_type + '\n'
            ind += 1
        out.append(line)
    with open(f"right_{ac1}",'w') as fp:
        fp.write(''.join(out))
    return f"right_{ac1}"

def update_prep_with_ac(ac_file, prep_file, output_file):
    """
    to correct the prep file's atom type from correct ac
    note that AC file is similar to pdb, it's column-width sensitive!
    """
    with open(prep_file, 'r') as f:
        prep_lines = f.readlines()
    for i in range(len(prep_lines)):
        if 'OMIT DU   BEG' in prep_lines[i]:
            res_name = prep_lines[i-1].split()[0]
            break
    n_atom, ac_mapping = 0, {}
    with open(ac_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith("ATOM"):
                n_atom += 1
                if line[17:21].strip() == res_name:
                    atom_name = line[12:17].strip()
                    atom_type = line.split()[-1]
                    ac_mapping[atom_name] = atom_type
    print(f"check ac_mapping and ac length: n_atom={n_atom}, ac_mapping={len(ac_mapping)}",flush=True)
    # prepi format: "%4d  %-5s %-5s %c%5d%4d%4d%10.3f%10.3f%10.3f%10.6f\n"
    # prepi LOOP: "%5s%5s\n"
    # prepi IMPROPER: "%5s%5s%5s%5s\n"
    # py: "{:4d}  {:<5s} {:<5s} {:1s}{:5d} ...
    updated_lines = []
    for line in prep_lines:
        parts = line.split()
        # get atom record lines
        if len(parts) >= 7 and parts[0].isdigit():
            atom_name = parts[1]
            prep_type = parts[2]
            if atom_name in ac_mapping:
                ac_type = ac_mapping[atom_name]
                if prep_type != ac_type:
                    _old = "{:<5s} {:<5s}".format(*parts[1:3])
                    parts[2] = ac_type
                    _new = "{:<5s} {:<5s}".format(*parts[1:3])
                    line = line.replace(_old, _new)
        updated_lines.append(line)

    with open(output_file, 'w', encoding='utf-8') as f:
        f.writelines(updated_lines)
    print(f"atom types updated in {output_file}",flush=True)

def get_lig_prep(prefix,resn,full_block_pdb,lig_indx):
    ## Risky! the lig is unsaturated in bonds! it will cause the error assignment of atom types!
    ## fix it by mapping the correct atom type from extend ac file!
    ## 20251231: extended ac file is still risky, because it just extended with one bond, use tmp/ac!
    cmds = [
    f"$AMBERHOME/bin/antechamber -fi pdb -i {prefix}.pdb -c rc -cf {prefix}.chg -fo ac -o {prefix}.ac -at gaff -seq no",
    f"$AMBERHOME/bin/prepgen -i {prefix}.ac -o tmp_{prefix}.prep -f int -rn {resn}",
    #f"$AMBERHOME/bin/parmchk2 -i {prefix}.mol2 -o {prefix}.frcmod -f mol2 -a Y"
    ]
    for cmd in cmds:
        print(f"running: {cmd}")
        os.system(cmd)
    right_ac = gen_right_ac(f"{prefix}.ac",'tmp/full_capped_resp_gaff.ac',lig_indx)
    update_prep_with_ac(right_ac, f"tmp_{prefix}.prep", f"{prefix}.prep")
    return f"{prefix}.prep"

def get_res_lib(prefix,resn,indices,FF="ff14SB"):
    for atom_type in ['gaff','amber']:
        cmd = f"$AMBERHOME/bin/antechamber -fi pdb -i {prefix}.pdb -c rc -cf {prefix}.chg \
               -fo mol2 -o {prefix}_{atom_type}.mol2 -at {atom_type}"
        print(f"running: {cmd}")
        os.system(cmd)
    ambMol2 = repalce_mol2_DU(f"{prefix}_amber.mol2",f"./tmp/full_capped_resp_amber.ac",indices)
    with open(f"{prefix}_amb.mol2",'w') as fp:
        fp.write(''.join(ambMol2))
    with open(f"{prefix}.tleap",'w') as fp:
        fp.write(f'''
        source leaprc.protein.{FF}
        {resn} = loadmol2 {prefix}_amb.mol2
        list
        desc {resn}
        #set {resn} head {resn}.{resn}.N
        #set {resn} tail {resn}.{resn}.C
        set {resn} head {resn}.1.N
        set {resn} tail {resn}.1.C
        set {resn} restype protein
        desc {resn}
        saveoff {resn} {prefix}.lib
        quit
        ''')
    os.system(f"$AMBERHOME/bin/tleap -f {prefix}.tleap")
    return f"{prefix}.lib"

def hybrid_mol2_types(amb_mol2,gaff_mol2,indicies):
    '''The two mol2 come from antechamber output using the same esp file with -seq no'''
    ind, flag, out = 0, False, []
    out_prefix = "full_capped_resp_hybrid"
    for aLine, gLine in zip(open(amb_mol2,'r'), open(gaff_mol2,'r')):
        if 'TRIPOS>ATOM' in aLine:
            flag = True
            out.append(aLine)
        elif 'TRIPOS>BOND' in aLine: flag = False
        if flag and 'TRIPOS>ATOM' not in aLine:
            if ind in indicies:
                out.append(gLine)
            else:
                out.append(aLine)
            ind += 1
        if not flag:
            out.append(aLine)
    with open(f"{out_prefix}.mol2",'w') as fp:
        fp.write(''.join(out))
    return out_prefix

def get_block_frcmod(pdb,indicies):
    '''getting correction connection from tmp/full_capped_resp_amber/gaff.ac, which is from esp, with
    optimized geometry!
    though in tmp dir, all residue names are fused (1 MOL residue), but that does not matter for frcmod file!
    parmchk2 -help:
    -s    ff parm set, it is suppressed by "-p" option
                      1 or gaff:    gaff (the default)
                      2 or gaff2:   gaff2
                      3 or parm99:  parm99
                      4 or parm10:  parm10
                      5 or lipid14: lipid14
                -frc  frcmod files to be loaded, the supported frcmods include
                      ff99SB, ff14SB, ff03 for proteins , bsc1, ol15, ol3 for DNA and yil for RNA
                      eg. ff14SB+bsc1+yil, ff99SB+bsc1
    '''
    for atp in ['amber','gaff']:
        pref = f'full_capped_resp_{atp}'
        #cmd0 = f"antechamber -fi ac -i tmp/{pref}.ac -fo mol2 -o {pref}.mol2 -c rc -cf resp2.chg -at {atp} -seq no"
        cmd0 = f"antechamber -fi ac -i tmp/{pref}.ac -fo mol2 -o {pref}.mol2 -c rc -cf resp2.chg -j 0 -seq no"
        print(f"running: {cmd0}",flush=True)
        res0 = subprocess.run(cmd0.split(),capture_output=True,text=True)
        print(f"stdout: {res0.stdout}\nstderr: {res0.stderr}\nreturncode: {res0.returncode}")

    if len(indicies) == 0 and type(indicies) == list:
        final_mol2 = "full_capped_resp_amber"
        param_sets = '1' #for lig-like residues with C/N/O, 3 can be wrong, e.g. getting a non-planar triazole moiety
    else:
        param_sets = '1'
        print(f"need to mix amber types with gaff types for {pdb}")
        final_mol2 = hybrid_mol2_types("full_capped_resp_amber.mol2",
                                       "full_capped_resp_gaff.mol2", indicies)
    cmd0 = f"parmchk2 -f mol2 -i {final_mol2}.mol2 -s {param_sets} -frc ff14SB -o {pdb}.frcmod -a Y"
    # maybe using resp2.chg and pymol exported mol2 is better for all blocks!
    print(f"running: {cmd0}")
    res = subprocess.run(cmd0.split(), capture_output=True, text=True)
    if res.returncode != 0:
        print(f"handling errors in parmchk2 {res.stderr}")
        print(f"errors occurred for mapping parameters, please check mol2 files in {os.getcwd()}!")
    #os.system(cmd1)
    return f"{pdb}.frcmod"

def get_bonding_by_id(leap_info,bond_info):
    desc_info, new_bond = {}, []
    flag = 0
    for line in leap_info.split('\n'):
        if line.startswith("Contents"):
            flag = 1
            continue
        elif "Quit" in line:
            break
        if flag == 1 and line.split()[0][2:] not in desc_info.keys():
            desc_info[line.split()[0][2:]] = []
        if flag == 1:
            desc_info[line.split()[0][2:]].append(line.split()[1].replace('>',''))
    capping_resi = []
    for x in ['ACE','NME']:
        if x in desc_info.keys():
            for y in desc_info[x]:
                capping_resi.append(y)
    for f in bond_info:
        # assuming the bonding info residues have all been unified in resname!
        i, j = f.split('.')[1], f.split('.')[3]
        new_bond.append(f.replace(i,desc_info[i][0]).replace(j,desc_info[j][0]))
    #removal = [f"remove cpx cpx.{y}\n" for y in desc_info['ACE']+desc_info['NME']]
    removal = [f"remove cpx cpx.{y}\n" for y in capping_resi]
    return new_bond, removal

def mapping_atom_back(xyz,pdb):
    # mapping the pdb atom name back to TS xyz
    # for PDB fed to antechamber, remove CONECT records
    lines = [x for x in open(xyz,'r')]
    atom_num = int(lines[0].strip('\n'))
    comment = lines[1]
    xyzs = lines[2:]
    i = 0
    news = [f"REMAKR {atom_num} atoms\nREMARK {comment}"]
    line2 = open(pdb,'r').readlines()
    atom_list = ['ATOM  ','HETATM']
    if len([x for x in line2 if x[0:6] in atom_list]) == 0:
        print(f"FATAL: no atoms in {pdb}")
        return None
    for line in line2:
        #if line.startswith("TER"): continue
        #if line.startswith('CONECT'): continue
        if line[0:6] not in atom_list:
            news.append(line)
        elif xyzs[i].split()[0] in line.strip('\n').split()[-1]:
            newxyz = [float(x) for x in xyzs[i].split()[1:]]
            newxyz = "{:8.3f}{:8.3f}{:8.3f}".format(*newxyz)
            oldxyz = line[30:54]
            #if line[17:20].strip() not in ['ACE', 'NME']:
            news.append(line.replace(oldxyz,newxyz))
            i += 1
        else:
            print(f"FATAL, atom orders aren not identical for"
                  f" {str(xyz)} {str(pdb)}, existing...")
            return None
    return news
    #with open(new,'w') as fp:
    #    fp.write(''.join(news))

def main():
    #tleap_files = []
    #pdb, esp, xyz = sys.argv[1], sys.argv[2], sys.argv[3]
    parser = argparse.ArgumentParser(
        description="post chg fitting procedures after building block slicing and optimization"
    )
    parser.add_argument('-pdb', '--pdb', dest='pdb',
                        help='pdb file with residue and atom information (Required)\n')

    parser.add_argument('-chg_type', '--charge_type', dest='chg_type', default='resp',
                        help='Charge method, either resp or bcc    (Required)'
                             f'e.g. -chg_type resp/bcc )\n')

    parser.add_argument('-xyz', '--xyz_file', dest='xyz',
                        help='the xyz file extracted from gaussian optimization\n')

    parser.add_argument('-sdf', '--sdf_file', dest='sdf',
                        help='the sdf file of the optimized block, required for bcc method\n')

    parser.add_argument('-esp', '--esp_file', dest='esp',
                        help=f'esp file from psi4/gaussian output, required for resp method\n')
    parser.add_argument('--capping_method', dest='capping_method',default='standard',
                        help=f'capping types from prepare_pep_caps, default = "standard"\n')
    parser.add_argument('--resname', dest='resname',
                        help=f"residue name must be given is capping_method is not standard\n")
    parser.add_argument('--charge_file',default='blocks_charges.json')
    args = parser.parse_args()
    base_pdb = '.'.join(args.pdb.split('.')[0:-1])
    opted_pdb = f"{base_pdb}_opt" if '_opt.' not in args.pdb else base_pdb
    if os.path.isfile('sqm.xyz'):
        args.xyz = 'sqm.xyz'
    if not os.path.isfile(f"{opted_pdb}.pdb"):
        with open(f"{opted_pdb}.pdb", 'w') as fp:
            mapped_pdb = mapping_atom_back(args.xyz, args.pdb)
            fp.write(''.join(mapped_pdb))
    # pdb_info : [resp0,ressInfo,ressType,bondinfo]
    pdb_info = get_init_chg(pdb=f"{opted_pdb}.pdb")
    init_chg = pdb_info[0]
    if args.chg_type == 'resp':
        all_resp = get_init_all_resp(args.esp,init_chg)
    elif args.chg_type == 'bcc':
        res_chgs = json.load(open(args.charge_file,'r'),encoding='utf-8')
        resn_chg = [res_chgs[x] for x in res_chgs.keys() if f'_{args.resname}.pdb' in x][0]
        all_resp = get_init_all_bcc(args.pdb,args.sdf,args.resname,resn_chg,pdb_info[4])
    ressInfo, ressType = pdb_info[1], pdb_info[2]
    ressList, liggList = [], []
    lig_indx = []
    for res in ressInfo.keys():
        indices = [x[0] for x in ressInfo[res]]
        charges = [all_resp[index] for index in indices]
        if ressType[res] == 'gaff':
            cur_lig_indx = [x[0] for x in ressInfo[res]]
            lig_indx += cur_lig_indx
            write_init_chg(charges,f"LIG_{res}.chg")
            with open(f'LIG_{res}.pdb', 'w') as ligF:
                ligF.write(''.join([x[1] for x in ressInfo[res]]))
            #TODO using indices is better also in ligand prep files!
            prepic = get_lig_prep(f'LIG_{res}',res,f"{opted_pdb}.pdb",cur_lig_indx)
            liggList.append(prepic)
        else:
            write_init_chg(charges, f"RES_{res}.chg")
            with open(f'RES_{res}.pdb', 'w') as resF:
                resF.write(''.join([x[1] for x in ressInfo[res]]))
            lib_file_name = get_res_lib(f'RES_{res}',res,indices)
            ressList.append(lib_file_name)

    frcmod = get_block_frcmod(pdb=opted_pdb,indicies=lig_indx)
    os.system("/usr/bin/find *.chk |/usr/bin/xargs rm -f")
    print("\n")
    print(f"tleap files:\n{' '.join(ressList)} {' '.join(liggList)} {frcmod}")

if __name__ == "__main__":
    main()
