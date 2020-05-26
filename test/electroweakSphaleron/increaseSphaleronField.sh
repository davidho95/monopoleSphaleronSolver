#PBS -lselect=2:ncpus=32:mem=2gb
#PBS -lwalltime=24:00:00

export LD_LIBRARY_PATH=/apps/gcc/8.2.0/lib64:$LD_LIBRARY_PATH
v=0.5
g=0.5
l=0.03125
q=0.286
b=1.2
module load mpi
inputDir="/rds/general/user/dlh16/home/monopoleSphaleronSolver/sphaleronDataPhysicalG0_5V0_5/sphaleronDataB0"
outputDir="/rds/general/user/dlh16/home/monopoleSphaleronSolver/sphaleronDataPhysicalG0_5V0_5"

mpiexec -n 64 /rds/general/user/dlh16/home/monopoleSphaleronSolver/bin/increaseSphaleronField -n 8 -m 8 -s 64 -p $outputDir -i $inputDir -g $g -v $v -l $l -q $q -b $b -B 1 -S 0 -I 5