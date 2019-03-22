#PBS -lselect=4:ncpus=32:mem=1gb
#PBS -lwalltime=24:00:00

export LD_LIBRARY_PATH=/apps/gcc/8.2.0/lib64:$LD_LIBRARY_PATH
module load mpi
mkdir $PBS_O_WORKDIR/fieldTestOutput
mpiexec -n 128 /rds/general/user/dlh16/home/monopoleSphaleronSolver/bin/fieldTest -s 128 -p $PBS_O_WORKDIR/fieldTestOutput -n 8 -m 16
cp -r $PBS_O_WORKDIR/fieldTestOutput /rds/general/user/dlh16/home/monopoleSphaleronSolver/output