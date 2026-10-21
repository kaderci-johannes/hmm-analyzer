if [ "$1" == "w" ]
then
watch -n 1 "condor_q okoseyan -nobatch;echo 'Output File:';ls -vah HTCondor/outputs/MC | grep '\.root';"
elif [ "$1" == "s" ]
then
watch -n 1 "condor_q okoseyan -nobatch | awk '{print \$1, \$2, \$5, \$6, \$11, \$NF}' | grep 'okoseyan'"
else
watch -w -n 1 "condor_q okoseyan -nobatch;echo 'Output File:';ls -vah HTCondor/outputs/MC | grep '\.root';"
fi
