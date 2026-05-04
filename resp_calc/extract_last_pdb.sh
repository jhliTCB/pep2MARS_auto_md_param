#!/bin/bash
input=$1
resName=$2
#abc=$(sed -n "/Stationary/, /Population/"p $input | grep -A 100 Coordinates | sed -n "/\-\-/, /\-\-\-/"p | sed /\-\-\-/d)
line=$(grep -n Coordinates $input | tail -1 | cut -f1 -d":")
#abc=$(sed -n ${line},$[line+80]p $input | sed '/[a-z]/d; /[A-Z]/d; /\-\-\-/d; /\=\=\=/d')
abc=$(sed -n "$[line+3],/\-\-\-\-\-/"p $input | sed '/[a-z]/d; /[A-Z]/d; /\-\-\-/d; /\=\=\=/d')
xyz=$(echo -e "$abc" | awk '{for (i=1; i<=NR; i++); \
if ($2=="1") $2="H"; \
if ($2=="6") $2="C"; \
if ($2=="7") $2="N"; \
if ($2=="8") $2="O"; \
if ($2=="15") $2="P"; \
if ($2=="16") $2="S"; \
if ($2=="17") $2="Cl"; \
if ($2=="26") $2="Fe"; \
if ($2=="9") $2="F"}; \
{printf " %-2s           %12.6f%12.6f%12.6f\n", $2,$4,$5,$6}')

echo -e "$xyz" > temp_file
line=$(cat temp_file | wc -l)
for i in $(seq 1 $line)
do
   atom_name=$(sed -n ${i}p temp_file | awk '{print $1}')
   XYZ=$(sed -n ${i}p temp_file | awk '{print $2,$3,$4}')
   if [[ $atom_name == 'Fe' ]]
   then
      FE='FE'
      echo -e "HETATM ${i} $FE $resName 1 $XYZ 0 0 $atom_name" |\
      awk '{printf "%-6s%5d %-5s%-4s%5d%12.3f%8.3f%8.3f%6.2f%6.2f%12s\n", $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}'
   else
      echo -e "HETATM ${i} $atom_name $resName 1 $XYZ 0 0 $atom_name" |\
      awk '{printf "%-6s%5d  %-4s%-4s%5d%12.3f%8.3f%8.3f%6.2f%6.2f%12s\n", $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}'
   fi
done


rm -f temp_file
