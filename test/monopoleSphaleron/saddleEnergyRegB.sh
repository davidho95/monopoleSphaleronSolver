#PBS -lselect=2:ncpus=32:mem=2gb
#PBS -lwalltime=24:00:00

export LD_LIBRARY_PATH=/apps/gcc/8.2.0/lib64:$LD_LIBRARY_PATH
module load mpi
mkdir /rds/general/user/dlh16/home/monopoleSphaleronSolver/saddle64Z2/saddleDataB12/saddleDataRegB
mpiexec -n 64 /rds/general/user/dlh16/home/monopoleSphaleronSolver/bin/saddleEnergyRegB -n 8 -m 8 -s 64 -i $PBS_O_WORKDIR/saddle64Z2/saddleDataB12/saddleData0_8 -p $PBS_O_WORKDIR/saddle64Z2/saddleDataB12/saddleDataRegB -g 0.5 -v 0.85 -l 0.125 -S 0.1 -E 0.66 -I 28 -b 1.5 -x 1
