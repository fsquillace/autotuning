

ID=$(qsub $1)

stat=$(qstat | grep feel | grep $ID | awk '{print $5}' )

while [ $stat != "C"  ];
do

  stat=$(qstat | grep feel | grep $ID | awk '{print $5}' )
  sleep 1

done

cat ex2.out | grep -e "^Time (sec):" -e "^Flops[:/]"
