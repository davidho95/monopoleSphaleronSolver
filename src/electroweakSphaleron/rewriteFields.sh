#PBS -lselect=2:ncpus=32:mem=2gb
#PBS -lwalltime=24:00:00

export LD_LIBRARY_PATH=/apps/gcc/8.2.0/lib64:$LD_LIBRARY_PATH
v=0.5
g=0.5
l=0.03125
q=0.286
module load mpi
dir1="/rds/general/user/dlh16/home/monopoleSphaleronSolver/sphaleronData/sphaleronDataPhysicalG0_5V0_5/sphaleronDataB0"
dir2="/rds/general/user/dlh16/home/monopoleSphaleronSolver/sphaleronData/sphaleronDataPhysicalG0_5V0_5/sphaleronDataB2"
dir3="/rds/general/user/dlh16/home/monopoleSphaleronSolver/sphaleronData/sphaleronDataPhysicalG0_5V0_5/sphaleronDataB4"

mpiexec -n 64 /rds/general/user/dlh16/home/monopoleSphaleronSolver/bin/rewriteFields -n 8 -m 8 -s 64 -p $dir1 -i $dir1 -g $g -v $v -l $l -q $q
mpiexec -n 64 /rds/general/user/dlh16/home/monopoleSphaleronSolver/bin/rewriteFields -n 8 -m 8 -s 64 -p $dir2 -i $dir2 -g $g -v $v -l $l -q $q
mpiexec -n 64 /rds/general/user/dlh16/home/monopoleSphaleronSolver/bin/rewriteFields -n 8 -m 8 -s 64 -p $dir3 -i $dir3 -g $g -v $v -l $l -q $q