#PBS -lselect=2:ncpus=32:mem=2gb
#PBS -lwalltime=24:00:00

export LD_LIBRARY_PATH=/apps/gcc/8.2.0/lib64:$LD_LIBRARY_PATH
B=0
v=1
g=1
l=0.125
q=0.1
module load mpi
outputDir="/rds/general/user/dlh16/home/monopoleSphaleronSolver/sphaleronData/sphaleronDataZ1Q0_1"
mkdir $outputDir

mpiexec -n 64 /rds/general/user/dlh16/home/monopoleSphaleronSolver/bin/generateSphaleron -n 8 -m 8 -s 64 -p $outputDir -g $g -v $v -l $l -q $q
