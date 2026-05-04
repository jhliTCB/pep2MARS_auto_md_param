#!/bin/bash
#final stage resp fitting
USAGE() {
    printf "a wrapper to fit resp charges from Gaussina/psi4's ESP
    USAGE:
    -pre|--prefix         the prefix name of output (.prep and .frcmod)
    -e|--esp              the .esp file name
    -i|--pdb              the template PDB file used for esp calculation
    -rn|--res-name        the residue name in .prep file (default=LIG)
    -pf|--purge-file      to clean up the tmp files from antechamber
    -pse|--pymol          provide a pse file for appending resp charge comparison (default='')
    -seq                  keep the input atom orders (without -seq will use the default tree orders
    -at|--atom-type       the forcefield atom types, default=gaff
    -j|--judge            bond order judge level in antechamber, see it's help for details
    -dr|--ac-doctor       with this option, it will turn off the ac doctor, useful for extacting .chg file from ac 
    \n"
    exit 0
}

get_prmtop() {
  cat > leap.in <<_EOF
source leaprc.gaff
loadamberparams ${1}.frcmod
loadamberprep   ${1}.prepc
saveamberparm $2 ${1}.prmtop ${1}.inpcrd
quit
_EOF
tleap -f leap.in
}

Jhelp() { printf "please specify 1-5 for atom/bond type judgement\n"; exit 1; }

UNKNOW() { printf "$1 is unkown, existing...\n"; exit 1; }
cleanup() {
rm -vf qout QOUT esout punch PREP.INF NEWPDB.PDB  *.AC *.AC0\
ATOMTYPE.INF ${prefix}*.ac leap.* #${prefix}*.chg
}

resN=LIG
purge=0
keep='' #-seq   atomic sequence order changable: yes(y)[default] or no(n)
##-at    atom type
##       gaff : the default
##       gaff2: for gaff2 (beta-version)
##       amber: for PARM94/99/99SB
##       bcc  : bcc 
##       sybyl: sybyl
atomtype='gaff'
custom_pse_file=''
glog=''
re='^[1-5]$'

while [[ $# -gt 0 ]]; do
  case $1 in
  -pre|--prefix)     prefix=$2;  shift 2;;
  -e|--esp)             esp=$2;  shift 2;;
  -l|--log)            glog=$2;  shift 2;;
  -i|--pdb)             pdb=$2;  shift 2;;
  -h|--help)            USAGE ;  shift 1;;
  -rn|--res-name)       resN=$2; shift 2;;
  -pf|--purge-file)     purge=1; shift 1;;
  -pse|--pymol)         custom_pse_file=$2; shift 2;;
  -seq)                 keep='-seq n'; shift 1;;
  -at|--atom-type)      atomtype=$2; shift 2;;
  -j|--judge)   [[ -z $2 ]] && Jhelp 
          [[ ! $2 =~ $re ]] && Jhelp; 
          Judge="-j $2";           shift 2;;
  -dr|--ac-doctor)      AcDoctor='-dr no'; shift 1;;
  *)                    UNKNOW $1       ;;
  esac
done
zz() { [[ -z $1 ]] && echo "no \$1 provided!" && exit 1; }
ZZ() { [[ $1 -ne 0 ]] && echo "something make the bash return non-zero please check" ; }
zz $prefix
[[ -z $esp && -z $pdb && $purge -eq 1 ]] && cleanup && exit 0
zz $esp
zz $pdb

#ml purge
#ml amber/a22t23

#sed -i "/CONECT/d" $prefix.pdb
cmds=(
"antechamber $Judge  $keep -fi gesp -fo ac -i $esp -o ${prefix}_resp.ac -c resp -ge ${prefix}.esp -pf y"
"antechamber $Judge $AcDoctor $keep -fi ac -i ${prefix}_resp.ac -c wc -cf ${prefix}.chg"
"[[ -f ${prefix}.chg ]] && echo yes"
"antechamber $Judge  $keep -fi pdb -i $pdb -c rc -cf ${prefix}.chg -fo ac -o ${prefix}_resp_pdb.ac -pf y"
"antechamber $Judge  $keep -fi pdb -i $pdb -c rc -cf ${prefix}.chg -fo mol2 -o ${prefix}_resp.mol2 -at sybyl -pf y"
"prepgen -i ${prefix}_resp_pdb.ac -o ${prefix}.prepc -f car -rn $resN"
"parmchk2 -i ${prefix}.prepc -o ${prefix}.frcmod -f prepc"
"get_prmtop ${prefix} $resN" )
for i in $(seq 1 ${#cmds[@]}); do
    ((j=i-1))
    echo "running cmd: ${cmds[$j]}"
    eval ${cmds[$j]}
    ZZ $?
done

if [[ ! -z $custom_pse_file ]]; then
    module load pymol/3.1.5.1_schrd
    cat > .pml <<_EOF
    cmd.load("${prefix}_resp.mol2")
    cmd.label("${prefix}_resp","'%.3f'%partial_charge")
    cmd.save("resp_$custom_pse_file")
    quit
_EOF
    $(which pymol) -c $custom_pse_file .pml
fi
[[ $purge -eq 1 ]] && cleanup
echo 0
