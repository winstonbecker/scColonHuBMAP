for value in {1,2,3,4,5,6,7,8,9}
do
for value2 in {1,2,3,4,5,6,7,8}
do
	echo $value2 
	echo $value
	sbatch 2b_runDecontX.sbatch $value2 $value
done
done

for value2 in {1,2,3,4}
do
	echo $value2 
	echo 10
	sbatch 2b_runDecontX.sbatch $value2 10
done