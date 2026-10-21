#!/bin/bash

for i in $@
do
	type=`echo $i | head -c2`
	name=`whoami`
	indir="/afs/cern.ch/work/${name:0:1}/$name/HTCondor/outputs"
	outdir="/afs/cern.ch/user/${name:0:1}/$name/projects/hmumu/test-analysis/build/src/drivers/HTCondor/outputs/Sums"
	folder="MC"

	if [ "$type" == "da" ]
	then
		folder="Data"
	fi

	if [ "$type" == "ef" ]
	then
		folder="Eff"
		outdir="."
		effs=($(ls -v ./HTCondor/outputs/Eff | grep $i | sed 's/_[0-9]*.root//'))
		eff0=${effs[0]}
		for eff in ${effs[@]}
		do
			echo $eff $eff0
			if [ "$eff" != "$eff0" ]
			then
				echo "Error: Mixed MCs"
				exit
			fi
		done
	fi

#	i=$eff0
	
	files=`ls -v $indir/$folder/$i*`

	if [ "$files" == "" ]
	then
		echo "No files found."
		exit
	fi
	
	echo
	ls -v $indir/$folder | grep $i
	echo -e "\nDo you wish to combine these files? (y/n)"
	
	while [ true ]
	read input
	do
		if [ "$input" == "y" ]
		then
			echo
			hadd -f $outdir/$i.root $files
			echo
			break
		elif [ "$input" == "n" ]
		then
			echo -e "\nNo file was created/changed.\n"
			break
		fi
	done
done
