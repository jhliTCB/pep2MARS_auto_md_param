#!/bin/bash -l
## The bash script here simply handle the the cases that having
## backbone-backbone cyclisation peptides, by copying the termini
## AAs to a new unit with differed names to those in Amber's 
## atuo termini residue mapping procedure, so that the cyclisation
## AAs in the last row of residues won't be treated as termini residues

tleap0=$1  # 'check_cpx.in'
complex=$2 # 'p4a.pdb'
tleap1=$3  # 'tleap_PP.in'

tleap -f $tleap0 >& leap0.log
#[[ $? -eq 0 ]] && exit 0
[[ $(grep 'at position' leap0.log|wc -l) -eq 0 ]] && exit 0

# 8.3f may make coordinates not matching to $complex!
positions=$(grep 'at position' leap0.log|sed 's/,//g; s/\.$//g'\
  |awk '{printf "%8.3f%8.3f%8.3f\n", $(NF-2),$(NF-1),$NF}'|sort|uniq)
line_num=$(echo -e "$positions"|wc -l)
cpx_lines=''
com_lines=''
for i in $(seq 1 $line_num); do
  coord=$(echo -e "$positions"|sed -n ${i}p)
  [[ $coord == '' || $coord == '\n' ]] && continue
  match0=$(grep "$coord" cpx_vac.pdb)
  match1=$(grep "$coord" $complex)
  cpx_lines="${cpx_lines}\n${match0}"
  com_lines="${com_lines}\n${match1}"
done
echo "position atoms in cpx_vac.pdb"
echo -e "$cpx_lines"
echo "position atoms in $complex"
echo -e "$com_lines"
cyclics=$(echo -e "$com_lines"|cut -c 18-26|sort|uniq)
cyc_len=$(echo -e "$cyclics"|wc -l)
to_be_copied=''
comp_new=cyclic_${complex}
cp -v  $complex $comp_new
for j in $(seq 1 $cyc_len); do
  oldResBlock=$(echo -e "$cyclics"|sed -n ${j}p)
  oldRes=$(echo -e "$oldResBlock"|awk '{print $1}')
  to_be_copied="${to_be_copied}\n${oldRes}"
  [[ $to_be_copied == '' || $to_be_copied == '\n' ]] && continue
  newResBlock=$(echo -e "$oldResBlock"|sed "s/$oldRes/${oldRes,,}/g")
  sed -i "s/ $oldResBlock/ ${newResBlock}/g" $comp_new
done

for res in $(echo -e "$to_be_copied"|sort|uniq); do
  [[ $res == '' || $res == '\n' ]] && continue
  #sed -i "/leaprc.protein/a  ${res,,} = copy ${res}" $tleap1
  # sometime the backbone linking is NCAA, need to copy just before cpx
  sed -i "/cpx = /i ${res,,} = copy ${res}" $tleap1
done
sed -i "s/${complex%.*}/${comp_new%.*}/g" $tleap1
