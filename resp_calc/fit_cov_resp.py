#!/public/home/lijunhao/soft/conda_envs/p4env/bin/python
# A small script to fit the resp charge for covalent system
# by Junhao Li
# The input must be a prepared complex with all hydrogen atoms!
import sys
import os
import subprocess
import math

global covResDict, ACE, NME, AMBTYPE
#covResDict: special residues, deprotonated with new charges
covResDict = {'CYS':'CYI','GLU':'GUU','SER':'SRR','LYS':'LYY',
              'HIS':'HII','HIP':'HII','HIE':'HII','HID':'HII',
              'ASP':'ASS','TYR':'TYY','THR':'THH'}
# the atomname
ACE = {'CH3':-0.3662,'H1':0.1123,'H2':0.1123,'H3':0.1123,
       'C':0.5972,'O':-0.5679,'CA':-0.3662,'HA':0.1123,
       'HA1':0.1123,'HA2':0.1123, 'HA3':0.1123, '1HH3':0.1123,
       'H01':0.1123,'H02':0.1123, 'H03':0.1123, 'H04':0.1123,
       '2HH3':0.1123, '3HH3':0.1123}
NME = {'N':-0.4157,'H':0.2719,'C':-0.1490,'CA':-0.1490,
       'H1':0.0976,'H2':0.0976,'H3':0.0976,'HA':0.0976,
       'H01':0.0976,'H02':0.0976,'H03':0.0976,'H04':0.0976,
       'HA1':0.0976,'HA2':0.0976,'HA3':0.0976, 'CH3':-0.1490,
       '1HH3':0.0976, '2HH3':0.0976, '3HH3':0.0976}
AMBTYPE = { # for pymol h_add, and for atomtype mapping for new-charge res
    'CYI':{'N':'N','H':'H','C':'C','O':'O','CA':'CX','CB':'2C',
    'SG':'SH','HB1':'H1','HB2':'H1','HB3':'H1','HA':'H1',
    'H01':'H1','H02':'H1','H03':'H1'
    },
    'GUU':{'N':'N','H':'H','C':'C','O':'O','CA':'CX','CB':'2C',
          'HB1':'HC','HB2':'HC','HB3':'HC','HA':'H1','CG':'2C',
          'HG1':'HC','HG2':'HC','HG3':'HC','CD':'CO','OE1':'O2',
          'OE2':'O2'
    },
    'SRR':{'N':'N','H':'H','C':'C','O':'O','CA':'CX','CB':'2C',
    'OG':'OH','HB1':'H1','HB2':'H1','HB3':'H1','HA':'H1',
    'H01':'H1','H02':'H1','H03':'H1'
    },
    'LYY':{'N':'N','H':'H','CA':'CX','HA':'H1','CB':'C8',
    'HB2':'HC','HB3':'HC','CG':'C8','HG2':'HC','HG3':'HC',
    'CD':'C8','HD2':'HC','HD3':'HC','CE':'C8','HE2':'HP',
    'HE3':'HP','NZ':'N3','HZ1':'H','HZ2':'H','HZ3':'H',
    'C':'C','O':'O'
    },
    'HII':{'CA':'CX','CB':'CT','C':'C','CD2':'CV','CD2':'CW',
    'CE1':'CR','CG':'CC','HA':'H1','HB2':'HC','HB3':'HC',
    'HD1':'H','HD2':'H4','HE1':'H5','HE2':'H','H':'H',
    'ND1':'NA','ND1':'NB','NE2':'NA','NE2':'NB',
    'N':'N','O':'O',
    },
    'ASS':{'N':'N','H':'H','CA':'CX','HA':'H1','CB':'2C',
    'HB2':'HC','HB3':'HC','CG':'CO','OD1':'O2','OD2':'O2',
    'C':'C','O':'O'
    },
    'TYY':{'N':'N','H':'H','CA':'CX','HA':'H1','CB':'CT',
    'HB2':'HC','HB3':'HC','CG':'CA','CD1':'CA','HD1':'HA',
    'CE1':'CA','HE1':'HA','CZ':'C','OH':'OH','HH':'HO',
    'CE2':'CA','HE2':'HA','CD2':'CA','HD2':'HA','C':'C','O':'O'
    },
    'THH':{'N':'N','H':'H','CA':'CX','HA':'H1','CB':'3C',
    'HB':'H1','CG2':'CT','HG21':'HC','HG22':'HC','HG23':'HC',
    'OG1':'OH','HG1':'HO','C':'C','O':'O'
    }
}
#def replace_cap_mc(line,atmNam):
#    if atmNam in ['N','C']:
#        new = atmNam+'1'
#        return line[0:12]+line[12:16].replace(atmNam+' ',new)+line[16:]
#    else:
#        return line

def repalce_mol2_type(inp,resid,covRes):
    ante = "%7d %-8s %10.4lf %10.4lf %10.4lf %-6s %5d %-8s %9.6lf\n"
    flag = False
    out = []
    for line in open(inp,'r'):
        if 'TRIPOS>ATOM' in line:   flag = True
        elif 'TRIPOS>BOND' in line: flag = False
        if flag and 'ATOM' not in line:
            atm = line.strip().split()
            if atm[6] == str(resid):
                newline = ante % (int(atm[0]),atm[1],float(atm[2]),float(atm[3]),
                float(atm[4]), AMBTYPE[covRes][atm[1]], int(atm[6]), atm[7],float(atm[8]))
                out.append(newline)
            else: out.append(line)
        else:
            out.append(line)
    return out

def get_init_chg(pdb,ligResName):
    # return a list masking the positions of capped atoms
    # read the REMARK of cov bond manully added to the _opt.pdb by psi4_esp.py
    # and split the pdb into reactive residue(s) and ligand
    ind, out, ress, resi, resn = 0, [], {}, {}, []
    resi[ligResName] = []
    lig = open('ligonly.pdb','w')
    new = open('ligress.pdb','w')
#    cap = []
    n_cap = 0
    for line in open(pdb, 'r'):
        if 'bondinfo:' in line: bondinfo = line.split(':')[-1]
        if line[0:6] in ['ATOM  ', 'HETATM']:
            atmNam = line[12:16].strip()
            resNam = line[17:20].strip()
            resNum = line[22:26].strip()
            if resNam != 'ACE' and resNam != 'NME':
                new.write(line[0:17].replace('HETATM','ATOM  ')+'LIS'+line[20:])
            if resNam == 'ACE':
                out.append(ACE[atmNam])
                n_cap += 1
#                cap.append(replace_cap_mc(line,atmNam))
            elif resNam in ['NME','NMA']:
                out.append(NME[atmNam])
#                cap.append(replace_cap_mc(line,atmNam))
                n_cap += 1
            else:
                out.append(0.0000)
            if resNam in covResDict.keys():
                if covResDict[resNam] not in ress.keys():
                    ress[covResDict[resNam]] = []
                    resi[covResDict[resNam]] = []
                    resn.append(resNum)
                ress[covResDict[resNam]].append(line.replace(resNam,covResDict[resNam]))
                resi[covResDict[resNam]].append(ind)
            if resNam == ligResName:
                lig.write(line)
                resi[ligResName].append(ind)
            ind += 1
    new.close()
    lig.close()
    ires = 0
    # Note that the CYI is a residue now, its mol2 can be written
    for res in ress.keys():
        with open('res{}.pdb'.format(str(ires)),'w') as resF:
            resF.write(''.join(ress[res]))
            #resF.write(''.join(cap))
        ires += 1
    print("number of capped atoms are {}".format(n_cap))
    return [out,resi,ress.keys(),resn,bondinfo]

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
    # charges are taken from the amber force field directly
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
def main():
    if len(sys.argv) < 4:
        print("usage: {} mol.pdb mol.esp ligResName".format(sys.argv[0]))
        print('put mol.pdb, mol.esp into current directory')
        exit(1)
    # for debug
    # exec(open('/public/software/apps/module/init/python.py').read())
    # module('purge')
    # module('load amber/a22t23 gaussian/16.c02')
    # module('list')
    pdb, esp = sys.argv[1], sys.argv[2]
    ligResName = sys.argv[3]
    #pdb, esp = 'mol.pdb', 'mol.esp'
    base_pdb = '.'.join(pdb.split('.')[0:-1])
    base_esp = '.'.join(esp.split('.')[0:-1])
    # get_init_chg will also unify and split the original MAE exported pdb:
    # ligress.pdb, ligonly.pdb, res0.pdb, res1.pdb, ... 
    pdb_info = get_init_chg(pdb=pdb,ligResName=ligResName)
    init_chg = pdb_info[0]
    write_init_chg(init_chg,fout='cap.chg')
    print(f"init_chg info for 2-step resp: {init_chg}")
    print(sum(init_chg),len(init_chg))
    
    os.makedirs('tmp',exist_ok=True)
    cmdx = f"$AMBERHOME/bin/espgen -i {esp} -o {base_esp}.dat"
    print(cmdx)
    os.system(cmdx)
    os.chdir('tmp')
    cmd0 = f"antechamber -fi gesp -i ../{esp} -fo ac -o {base_esp}.ac -c resp"
    os.system(cmd0)
    os.chdir('..')
    
    mod_resp_inp(init_chg,'tmp/ANTECHAMBER_RESP1.IN',stage=1,fout='resp1.in')
    mod_resp_inp(init_chg,'tmp/ANTECHAMBER_RESP2.IN',stage=2,fout='resp2.in')
    cmd1 = f"$AMBERHOME/bin/resp -O -i resp1.in -o resp1.out -p resp1.pch -t resp1.chg -q cap.chg -e {base_esp}.dat"
    cmd2 = f"$AMBERHOME/bin/resp -O -i resp2.in -o resp2.out -p resp2.pch -t resp2.chg -q resp1.chg -e {base_esp}.dat"
    print(f"running custom resp fit\n{cmd1}\n{cmd2}\n")
    os.system(cmd1)
    os.system(cmd2)
    # now the resp2.chg can be used, the sum charges of the capped atoms should zero.
    
    #TODO fix the cases with more than one covalent bound residues
    # 
    deprot = [x for x in pdb_info[2] if x != ligResName]
    resi = pdb_info[1]
    ress = {'ligress':'LIS','ligonly':ligResName,'res0':deprot[0]}
    all_resp = []
    for line in open('resp2.chg','r'):
        all_resp += [float(x) for x in line.strip().split()]
    
    ligonly = [all_resp[x] for x in resi[ress['ligonly']]]
    res0    = [all_resp[x] for x in resi[ress['res0']]]
    covRes  = deprot[0] # to be fixed
    ligress = [all_resp[x] for x in resi[ress['res0']]+resi[ress['ligonly']]]
    print(f"ligand atoms: {len(ligonly)}\n{deprot[0]} atoms:  {len(res0)}\nTotal atoms:  {len(all_resp)}")
    print(f"ligand charge {sum(ligonly)}\n{deprot[0]} charge: {sum(res0)}\nall_resp charge: {sum(all_resp)}")
    print(f"ligress charge {sum(ligress)}\n{ligress} atoms: {len(ligress)}\n")
    write_init_chg(ligonly,fout='ligonly.chg')
    write_init_chg(res0, fout='res0.chg')
    write_init_chg(ligress, fout='ligress.chg')
    
    for prefix in ['ligress','ligonly','res0']:
        res = ress[prefix]
        atm, dr = 'gaff', '  '
        if 'lig' not in prefix: atm='amber'
        if prefix == 'ligress': dr='-dr no'
        cmd3 = f"$AMBERHOME/bin/antechamber -fi pdb -i {prefix}.pdb -c rc -cf {prefix}.chg -fo ac -o {prefix}.ac -at {atm} {dr}"
        cmd4 = f"$AMBERHOME/bin/antechamber -fi pdb -i {prefix}.pdb -c rc -cf {prefix}.chg -fo mol2 -o {prefix}.mol2 -at {atm} {dr}"
        cmd5 = f"$AMBERHOME/bin/prepgen -i {prefix}.ac -o {prefix}.prep -rn {res}"
        cmd6 = f"$AMBERHOME/bin/parmchk2 -i {prefix}.mol2 -o {prefix}.frcmod -f mol2 -a Y"
        print(f"running:\n{cmd3}\n{cmd4}\n{cmd5}\n{cmd6}\n")
        os.system(cmd3)
        os.system(cmd4)
        if prefix == 'ligonly': os.system(cmd5)
        else:
            newlines = repalce_mol2_type(f'{prefix}.mol2',pdb_info[3][0],covRes)
            with open(f'{prefix}.mol2','w') as newmol2:
                newmol2.write(''.join(newlines))
        if prefix == 'ligress':
            with open('tmp.in','w') as fp:
                fp.write(f'p = loadmol2 {prefix}.mol2\n')
                fp.write(f'bond {pdb_info[4]}\n')
                fp.write(f'savemol2 p {prefix}.mol2 1\nquit\n')
            os.system('tleap -f tmp.in')
        os.system(cmd6)
    
    FF='source leaprc.protein.ff14SB'
    res0=deprot[0]
    with open('res0.in','w') as leap1:
        leap1.write(f'''{FF}
    {res0} = loadmol2 res0.mol2
    list
    desc {res0}
    set {res0} head {res0}.{res0}.N
    set {res0} tail {res0}.{res0}.C
    set {res0} restype protein
    desc {res0}
    saveoff {res0} res0.lib
    quit''')
    
    os.system("$AMBERHOME/bin/tleap -f res0.in")
    # for better visualization:
    os.system('antechamber -fi mol2 -i ligress.mol2 -fo mol2 -o ligress_sybyl.mol2 -at sybyl -dr no')
    print("Files to be used for covMD: res0.lib, ligonly.prep ligress.frcmod")
    print("Files to be used for deprotonated CYS MD: res0.lib, ligonly.prep, ligonly.frcmod")
    
if __name__ == "__main__":
    main()