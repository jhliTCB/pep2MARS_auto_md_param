#!/bin/bash
# This is the top-level wrapper for the entire workflow

set -e

usage() {
  cat <<EOF
Usage:
  $0 \\
    -w|-workDir WORKDIR \\
    -p|-pdbName PDBNAME \\
    -j|-jobName JOBNAME \\
    -t|-jobTime JOBTIME \\
    -f|-forceField FORCEFIELD \\
    -T|-temper TEMPER \\
    -b|-waterBox WATERBOX \\
    -c|-chargeType CHARGETYPE \\
    -s|-saltConc SALTCONC \\
    [-d|-scriptDir SCRIPTDIR]

Example:
  $0 \\
    -w /path/to/work \\
    -p protein.pdb \\
    -j test1 \\
    -t 24 \\
    -f amber_ff14SB \\
    -T 300 \\
    -b TIP3PBOX \\
    -c bcc \\
    -s 0.15

Options:
  -w, -workDir       Working directory
  -p, -pdbName       PDB file name
  -j, -jobName       Job name
  -t, -jobTime       Job time
  -f, -forceField    Force field
  -T, -temper        Temperature
  -b, -waterBox      Water box
  -c, -chargeType    Charge type
  -s, -saltConc      Salt concentration
  -d, -scriptDir     Script directory
  -h, -help          Show this help message
EOF
}

scriptDir="${scriptDir:-$(dirname "$(realpath "$0")")}"
[[ -z $1 ]] && usage && exit 0
while [[ $# -gt 0 ]]; do
  case "$1" in
    -w|-workDir) workDir="$2"; shift 2 ;;
    -p|-pdbName) pdbName="$2"; shift 2 ;;
    -j|-jobName) jobName="$2"; shift 2 ;;
    -t|-jobTime) jobTime="$2"; shift 2 ;;
    -f|-forceField) forceField="$2"; shift 2 ;;
    -T|-temper) temper="$2"; shift 2 ;;
    -b|-waterBox) waterBox="$2"; shift 2 ;;
    -c|-chargeType) chargeType="$2"; shift 2 ;;
    -s|-saltConc) saltConc="$2"; shift 2 ;;
    -d|-scriptDir) scriptDir="$2"; shift 2 ;;
    -h|-help) usage; exit 0 ;;
    *) echo "Error: unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

missing=0

for var in workDir pdbName jobName jobTime forceField temper waterBox chargeType saltConc; do
  if [[ -z "${!var}" ]]; then
    echo "Error: missing required argument: $var" >&2
    missing=1
  fi
done

if [[ "$missing" -eq 1 ]]; then
  echo
  usage
  exit 1
fi

echo "arguments:"
echo "workDir=$workDir"
echo "pdbName=$pdbName"
echo "jobName=$jobName"
echo "jobTime=$jobTime"
echo "forceField=$forceField"
echo "temper=$temper"
echo "waterBox=$waterBox"
echo "chargeType=$chargeType"
echo "saltConc=$saltConc"
echo "scriptDir=$scriptDir"

ff=$(echo "$forceField" | cut -f2 -d'_')
echo "ff=$ff"
topDir=$PWD
mkdir -p $workDir
cd "$workDir"
echo "forceField=$forceField"

if [[ "$jobTime" =~ \. ]]; then
  hours=10 # 10 hours for parameter fitting
  minutes=$(echo "$jobTime" | awk '{printf "%d", $1*60}')
else
  hours="$jobTime"
  minutes="00"
fi

echo "gonna run $hours hours, $minutes minutes"

cat > master_worker.sh <<_EOF
#!/bin/bash
# Use only 1 CPUs here, because ESP charge will be Sbatch
# in build_system.py with 24 CPUs
#SBATCH --job-name="prep_$jobName"
#SBATCH -p CPU_5318Y_96C
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t $hours:$minutes:00

ml purge
ml amber/a22t23
set -e

echo Preparing the Amber prmtop and inpcrd files for $pdbName

$(which python) $scriptDir/1_build_system.py "$pdbName" "blank_string" \\
 "$chargeType" "$ff" "$temper" "$waterBox" "$jobTime" "$saltConc"

_EOF

sbatch master_worker.sh
cd $topDir
