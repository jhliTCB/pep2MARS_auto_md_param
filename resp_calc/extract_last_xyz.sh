#!/bin/bash
input=$1
chgmul=$(grep -i multiplicity $input|uniq|tail -1|awk '{print $3,$NF}')
#abc=$(sed -n "/Stationary/, /Population/"p $input | grep -A 100 Coordinates | sed -n "/\-\-/, /\-\-\-/"p | sed /\-\-\-/d)
line=$(grep -n Coordinates $input | tail -1 | cut -f1 -d":")
#abc=$(sed -n ${line},$[line+80]p $input | sed '/[a-z]/d; /[A-Z]/d; /\-\-\-/d; /\=\=\=/d')
abc=$(sed -n "$[line+3],/\-\-\-\-\-/"p $input | sed '/[a-z]/d; /[A-Z]/d; /\-\-\-/d; /\=\=\=/d')
xyz=$(echo -e "$abc" | awk '{for (i=1; i<=NR; i++); \
if ($2=="1") $2="H"; \
if ($2=="6") $2="C"; \
if ($2=="7") $2="N"; \
if ($2=="8") $2="O"; if ($2=="15") $2="P"; \
if ($2=="16") $2="S"; \
if ($2=="17") $2="Cl"; \
if ($2=="26") $2="Fe"; \
if ($2=="9") $2="F"; \
if ($2=="35") $2="Br"; \
if ($2=="53") $2="I"}; \
{printf " %-2s           %12.6f%12.6f%12.6f\n", $2,$4,$5,$6}')
get_coord() { echo -e "$xyz"|sed -n ${1}p|awk '{print $2,$3,$4}'; }
[[ ! -z $2 && ! -z $3 ]] && dist=$(measure.py $(get_coord $2) $(get_coord $3)) \
&& dist="dist:$2,$3:$dist" || dist=""
[[ ! -z $4 && ! -z $5 && ! -z $6 ]] && angl=$(measure.py $(get_coord $4) $(get_coord $5) \
$(get_coord $6)) && angl="angl:$4,$5,$6:$angl" || angl=""
ener_term="freq0,freq1,freq2,ZPE,ener:$(get_freq_energy.sh $input|grep "log\|out" | awk '{print $3,$4,$5,$6,$7}'|sed 's/ /,/g')"
num_atom=$(echo -e "$xyz"|wc -l)
echo "$num_atom"
echo "$input;chgmul:$chgmul;$dist;$angl;$ener_term"
echo -e "$xyz"
