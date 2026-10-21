#!/bin/bash
isNotInsideLoop_mc=true
isNotInsideLoop_data=true
processes=("$@")
R=${@: -1}
if [ "$R" == "R" ]
then
	unset 'processes[-1]'
fi

for i in "${processes[@]}"
do
	echo
	f_2=$(head -c2 <<< "$i")
	for f in `ls -v HTCondor/lists/ | grep "$i.*\.list" | sed 's/\.list//'`
	do
		if [ "$f_2" == "mc" ]
		then
			N=()
			eff=$(sed 's/_[0-9]*$//' <<< eff${f:2}) #"eff${i:2}"
			f_lxsw=$(sed 's/_[0-9]*$/.lxsw/' <<< "./HTCondor/lists/$f")
			for n in `cat $f_lxsw`
			do
				N+=($n)
			done
			echo "true,${N[0]},${N[1]},${N[2]},$eff,$f"
		elif [ "$f_2" == "da" ]
		then
			echo "false,"$f
		fi
	done
	
	if [ "$R" == "R" ]
	then
		echo -e "\nRelevant joblist(s) will be overwritten with the lines above."
	else
		echo -e "\nThe lines above will be appended to their relevant joblist(s)."
	fi
	echo "Are you sure? (y/n)"
	
	while [ true ]
	read input
	do
		if [ "$input" == "y" ]
		then
			for f in `ls -v HTCondor/lists/ | grep "$i.*\.list" | sed 's/\.list//'`
			do
				if [ "$f_2" == "mc" ]
				then
					N=()
					eff=$(sed 's/_[0-9]*$//' <<< eff${f:2}) #"eff${i:2}"
					f_lxsw=$(sed 's/_[0-9]*$/.lxsw/' <<< "./HTCondor/lists/$f")
					for n in `cat $f_lxsw`
					do
						N+=($n)
					done
					if [ "$R" == "R" ] && $isNotInsideLoop_mc
					then
						rm joblist_mc.txt
						rm joblist_eff.txt
						isNotInsideLoop_mc=false
					fi
					echo "true,${N[0]},${N[1]},${N[2]},$eff,$f" >> joblist_mc.txt
					echo "true,${N[0]},${N[1]},${N[2]},eff${f:2},$f" >> joblist_eff.txt
				elif [ "$f_2" == "da" ]
				then
					if [ "$R" == "R" ] && $isNotInsideLoop_data
					then
						rm joblist_data.txt
						isNotInsideLoop_data=false
					fi
					echo "false,"$f  >> joblist_data.txt
				fi
			done
			if [ "$R" == "R" ]
			then
				echo "Relevant joblist(s) were recreated."
				R=""
			else
				echo "List files were appended to their relevant joblist(s)"
			fi
			break
		elif [ "$input" == "n" ]
		then
			echo Joblists were not changed.
			break
		fi
	done
done
