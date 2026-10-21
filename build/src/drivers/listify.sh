#!/bin/bash
voms-proxy-info || voms-proxy-init -voms cms -valid 192:00

flag=""
if [ "${@: -1}" == "A" ]
then
	flag="-aggregate"
fi

F=()
datasets=$(sed 's/,/ /g' <<< "$1")
for dataset in $datasets
do
	F+=$(env -u PYTHONHOME -u PYTHONPATH dasgoclient -query "file dataset=$dataset" $flag)
	F+=$'\n'
done

i=0
j=1
for f in $F
do
	if [ $(($i%10)) -eq 0 ]
	then
		rm HTCondor/lists/$2_$j.list 2> /dev/null
	fi
	echo "root://cms-xrd-global.cern.ch/$f" >> HTCondor/lists/$2_$j.list
	((i++))
	if [ $(($i%10)) -eq 0 ]
	then
		((j++))
	fi
done
if [[ "$3" != "" && "$4" != "" ]]
then
	root -l -q 'LXSW.C("'$2'",'$3','$4')'
fi
