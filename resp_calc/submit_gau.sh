#!/bin/bash
# A thin gaussian sbatch wrapper, by Junhao Li
CDEF="\033[0m"     ;  b_CDEF="\033[1m"      # default; bold_default
CRED="\033[0;31m"  ;  b_CRED="\033[1;31m"   # red;     bold_red
CGRE="\033[0;32m"  ;  b_CGRE="\033[1;32m"   # green;   bold_green
CYEL="\033[0;33m"  ;  b_CYEL="\033[1;33m"   # yellow;  bold_yellow
CBLU="\033[0;34m"  ;  b_CBLU="\033[1;34m"   # blue;    bold_blue
CMAG="\033[0;35m"  ;  b_CMAG="\033[1;35m"   # magenta; bold_magenta
CCYA="\033[0;36m"  ;  b_CCYA="\033[1;36m"   # cyan;    bold_cyan

print_help() {
printf "# A thin wrapper to submit gaussian jobs, please ensure that the number of\n"
printf "# CPU and memory are specific in the .com/.gjf file, e.g.\n"
printf " ${b_CRED}%%nproc=16${CDEF}\n"
printf " ${b_CRED}%%mem=8GB${CDEF}\n\n"
printf "options:\n"
printf " ${b_CDEF}-h|--help${CDEF}       print this information\n"
printf " ${b_CDEF}-i|--input${CDEF}      provide the gaussian input file, ${b_CGRE}required${CDEF}\n"
printf " ${b_CDEF}-j|--jobname${CDEF}    provide the jobname for slurm, ${b_CMAG}default: input file prefix${CDEF}\n"
printf " ${b_CDEF}-c|--core${CDEF}       proivde the number of cores(threads) for slurm,\n" 
printf "                ${b_CMAG}default: reading from .com/.gjf${CDEF}\n"
printf " ${b_CDEF}-m|--memory${CDEF}     provide the memobry usage for slurm,\n"
printf "                ${b_CMAG}default: reading from .com/.gjf${CDEF}\n"
printf " ${b_CDEF}-t|--time${CDEF}       provide the time length of the slurm job, ${b_CMAG}default=80:90:00${CDEF}\n"
printf " ${b_CDEF}-d|-dry${CDEF}         just generates the slurm script (${b_CGRE}ttt.sh${CDEF}) without submission;\n"
printf "                 default is to submit it immediately\n"
printf " ${b_CDEF}-w|--node${CDEF}       provide the name of the computing node, e.g. ${b_CBLU}-w comput8,${CDEF}\n"
printf "                 default: not provide, let the slurm decide\n\n"
printf "\n To submit multiple .com/.gjf files, use a bash for clause, e.g.\n"
printf " ${b_CGRE}for f in *.com; do $1 -i \$f; done${CDEF}\n\n"
exit
}

give_input() { echo "The Gaussian input must be given!\nUSAGE:";}

sub() {
cat <<_EOF
#!/bin/bash
#SBATCH -J $1
#SBATCH -N 1
#SBATCH -c $2
#SBATCH --time $3
#SBATCH -p $6
#SBATCH --mem $4
$7

ml purge
ml gaussian/16.c02
g16 $5
[[ ! -z \$(ls \$GAUSS_SCRATCH|grep esp) ]] && cp -vf \$GAUSS_SCRATCH/*.esp . 2> /dev/null
exit 0
_EOF
}

INP=none
DRY=0

[[ $# == 0 ]] && print_help $0

while [[ $# -gt 0 ]]; do
    case $1 in
        -j | --jobname)        JOBN=$2   ;        shift 2 ;;
        -c | --core)           CORE=$2   ;        shift 2 ;;
        -t | --time)           TIME=$2   ;        shift 2 ;;
        -m | --mem)            MEMO=$2   ;        shift 2 ;;
        -i | --input)           INP=$2   ;        shift 2 ;;
        -p | --partition)      PART=$2   ;        shift 2 ;;
        -d | -dry | --dry_run)  DRY=1    ;        shift 1 ;; 
        -w | --node)           node=$2   ;        shift 2 ;;
        -h | --help)        print_help $0;        shift 1 ;;
        *)                  print_help $0;        shift 1 ;;
    esac
done

[[ $INP == 'none' ]] && give_input && print_help $0
[[ -z $JOBN ]] && JOBN=${INP%.*}
[[ -z $CORE ]] && CORE=$(grep -i "proc\|cpu" $INP|grep "%"|cut -f2 -d'='|uniq)
[[ -z $MEMO ]] && MEMO=$(grep -i "mem" $INP|grep "%"|cut -f2 -d'='|sort|tail -1)
[[ -z $TIME ]] && TIME="80:90:00"
[[ -z $CORE || -z $MEMO ]] && echo cpu and/or mem not speicified && print_help $0
[[ -z $PART ]] && PART=CPU_5318Y_96C
[[ -z $node ]] && NODE="###" || NODE="#SBATCH -w $node"

sub $JOBN $CORE $TIME $MEMO $INP $PART "$NODE" > ttt.sh
[[ $DRY == 0 ]] && sbatch ttt.sh
