#PBS -lselect=2:ncpus=32:mem=2gb
#PBS -lwalltime=24:00:00

export LD_LIBRARY_PATH=/apps/gcc/8.2.0/lib64:$LD_LIBRARY_PATH
module load mpi
mkdir $PBS_O_WORKDIR/saddle64G0_5/finiteSize48
mpiexec -n 64 /rds/general/user/dlh16/home/monopoleSphaleronSolver/bin/findSaddle -n 8 -m 8 -s 48 -p /rds/general/user/dlh16/home/monopoleSphaleronSolver/output/saddleDataG0_5B8/finiteSize48
