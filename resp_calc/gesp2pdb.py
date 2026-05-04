#!/usr/bin/env python3
# convert the charges and coordinates of a ".esp" file to pdb
# for visulization of the ESP charges and fitting points
# by Junhao Li
import sys
def readPDB(pdb):
    out, con = [], []
    for line in open(pdb,'r'):
        if line[0:6] in ['ATOM  ', 'HETATM']:
            out.append([line[0:6],line[6:11],line[12:16],line[16:17],
                line[17:20],line[21:22],line[22:26],line[26:27],
                line[30:38],line[38:46],line[46:54],line[54:60],
                line[60:66],line[76:78],line[78:80]])
        if 'CONECT' in line or 'END' in line: con.append(line)
    return [out, con]

def coords2pdb(data,tmpl='None'):
    PDBForm = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}" + \
            "{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n"
    outList = []
    b2a = 0.52918
    if tmpl != 'None': oriPDB = readPDB(tmpl)
    for i in range(len(data)):
        reclst = data[i].split()
        if len(reclst) == 5:
            elemen = reclst[0]
            resNam = 'LIG'
            chgInd = -1
            if tmpl == 'None':
                atmNam = reclst[0]+str(i+1)
                resNam = 'LIG'
            else:
                atmNam = oriPDB[0][i][2]
                resNam = oriPDB[0][i][4]
        else:
            elemen = "pts"
            resNam = "PTS"
            atmNam = "PTS"
            chgInd = 0
        coords = [b2a*float(x.replace('D','E')) for x in reclst[1:4]]
        resNum = 999
        atmNum = i+1
        charge = float(reclst[chgInd].replace('D','E'))
        atmlst = ['ATOM',i+1,atmNam,' ',resNam,"X",resNum," ",]+coords+\
                 [1.00,charge,elemen,"  "]
        outList.append(PDBForm.format(*atmlst))
    if tmpl != 'None': outList += oriPDB[1]
    return outList

def get_coords(gesp):
    lig, pts, rmk = [], [], []
    fl1, fl2 = False, False
    for line in open(gesp):
        if 'ATOMIC COORDINATES AND ESP CHARGES' in line:
            fl1 = True
            n_atoms = int(line.split()[-1])
        if fl1 and 'ATOMIC COORDINATES' not in line:
            if 'DIPOLE MOMENT' not in line:
                lig.append(line)
            else:   fl1 = False
        if 'ESP VALUES AND GRID' in line:
            fl2 = True
            n_points = int(line.split()[-1])
        if fl2 and 'ESP VALUES' not in line:
            if len(line) > 10:
                pts.append(line)
            else: fl2 = False
        if not fl1 and not fl2:
            rmk.append(f"REMARK {line}")
    return [rmk,lig,pts,n_atoms,n_points]

def writePDB(pdblst,out):
    outF = open(out, 'w')
    outF.write(''.join(pdblst))
    outF.close()

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("usage: {} orig.lig.pdb gau/psi4.esp ".format(sys.argv[0]))
        exit(1)
    name = sys.argv[2][:-4]
    data = get_coords(sys.argv[2])
    writePDB(data[0]+coords2pdb(data[1],tmpl=sys.argv[1]),f"{name}.lig.pdb")
    writePDB(data[0]+coords2pdb(data[2]),f"{name}.esp.pdb")
