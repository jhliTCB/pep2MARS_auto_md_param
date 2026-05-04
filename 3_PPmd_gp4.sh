#!/bin/bash
#cmd = [f"{scriptDir}/md/P_P_MD/gp4.sh", os.getcwd(),
#               "-wd",wd,"-pt","lig","-inp",f,"-opt","gau",
#               "-chg",chg,"-ct",ct,"--resn",resn, "--cpus",24, "--mems",xx]
topDir=$PWD
workDir=$1
jobName=$3
cd $workDir
echo "pdb file in $PWD: $(ls *.pdb)"
all_args=($@)
charge_fitting_args=${@:2}
while [[ $# -gt 0 ]]; do
  case $1 in
    -wd) chargeDir=$2; shift 2;;
    -inp) inp_file=$2; shift 2;;
    -chg)  chg_val=$2; shift 2;;
    --mems) mem_val=$2; shift 2;;
    --cpus) cpu_val=$2; shift 2;;
    *) shift 1;;
  esac
done
echo "charge fitting dir is : $chargeDir"
echo "args: ${all_args[@]}"
echo "shifted args: $@"
# Note that ${0%/*} requires a path, so, don't add it to $PATH
scriptDir=${0%/*}

if [[ $(hostname) =~ 192.168 ]]; then
  python_exec1='/Users/envs/amb22/bin/python'
else
  python_exec1='/public/home/lijunhao/soft/conda_envs/fldev/bin/python'
  python_exec2='/public/home/lijunhao/soft/conda_envs/p4dev/bin/python'
fi

new_inp=${inp_file%.*}_opt.pdb
mkdir -p $chargeDir
cp -vf ${inp_file%.*}.* $chargeDir/
just_opt=$(echo ${charge_fitting_args}|sed "s/ bcc / resp /g")
new_arg=$(echo ${charge_fitting_args}|sed \
   "s/$inp_file/$new_inp/g; s/opt gau/opt none/g; s/wd $chargeDir/wd \$PWD/g")
#new_arg=$(echo ${charge_fitting_args}|sed "s/opt gau/opt none/g")

cd $chargeDir
cat > chargeFitting.sh <<_EOF
#!/bin/bash
# another worker due to failures in openMP running using subprocessing in Ubuntu
#SBATCH -J $jobName
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem 800MB
#SBATCH -t 32:00:00
#SBATCH -p CPU_5318Y_96C

inp_file=$inp_file
inp_pref=${inp_file%.*}
opt_pref=\${inp_pref}_opt
original_args="${all_args[@]}"
new_arg="${new_arg}"

cmd="$python_exec1 $scriptDir/utils/mm_opt.py -ipdb \$inp_file -isdf \${inp_pref}.sdf \\
      -opdb \${opt_pref}.pdb -osdf \${opt_pref}.sdf  --maxIts 3000"
  echo Running simple mm opt:
  echo \$cmd
  eval "\$cmd"
  if [[ \$? -ne 0 ]]; then
    echo "first attempt of opt failed, tried without caps relax"
    cmd2="\${cmd} --no_pre_relax_caps"
    echo "\$cmd2"
  fi
  eval "\$cmd2"
  if [[ \$? -ne 0 ]]; then
    echo "second attempt of opt failed! tried without trans_amide constraints!"
    cmd3="\${cmd2} --no_trans_constraints"
    echo "\$cmd3"
  fi
  eval "\$cmd3"


  [[ \$inp_file =~ "pure_lig_" ]] || new_arg="\${new_arg} --ac_flags '@dr no'"

cat > .tmp.\${inp_pref}.sh <<_EOE
#!/bin/bash
#SBATCH -J $jobName
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c $cpu_val
#SBATCH --mem ${mem_val}GB
#SBATCH -t 32:00:00
#SBATCH -p CPU_5318Y_96C
###SBATCH -p GPU_4090_8
cmd="$python_exec2 $scriptDir/resp_calc/atomic_charge_fitting.py \$new_arg"
echo Running atomic_charge_fitting cmd:
echo \\\$cmd
eval "\\\$cmd"
_EOE

jid=\$(sbatch .tmp.\${inp_pref}.sh|awk '{print \$NF}')
sleep 5
echo jid for resp/bcc using _opt.pdb is \$jid
#sacct_out=\$(sacct -X --parsable2 --noheader --format JobID,State,Elapsed,End -j \$jid)
#while [[ ! -z \$(squeue -j \$jid 2> /dev/null) ]]; do
ii=0
while [[  ! -z \$(squeue -j \$jid 2> /dev/null| grep -v "JOBID PARTITION     NAME")  ]]; do
  sleep 60
  ((ii+=1))
  [[ \$ii -gt 720 ]] && echo 12h passed, breaking it anyway! && break
done


echo charge_fitting jobDir is $chargeDir
#cp -vf \${jobDir}/*.esp \${jobDir}/*.pdb \${jobDir}/*.dat \${jobDir}/*.lib . 2> /dev/null
#cp -vf \${jobDir}/*.pse \${jobDir}/*.prep* \${jobDir}/*.frcmod . 2> /dev/null

find -type f  | xargs chmod 660
find -type d | xargs chmod 770
echo 0
_EOF

sbatch chargeFitting.sh
cd $topDir
