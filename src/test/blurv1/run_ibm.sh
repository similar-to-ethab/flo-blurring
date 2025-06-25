#!/bin/bash
#SBATCH -J FBA_3D_icx		# Job name
#SBATCH -o fba3d.o%j		# Name of stdout output file
#SBATCH -e fba3d.e%j		# Name of stderr error file
#SBATCH -p skx-dev  		# Queue (partition) name
#SBATCH -N 1			dd# Total # of nodes
#SBATCH -n 1  		# Total # of MPI ranks (cores)
#SBATCH -t 48:00:00		# Run time (hh:mm:ss)
#SBATCH --mail-type=all		# Send email at begin and end of job
#SBATCH --mail-user=ethanps@email.sc.edu
#SBATCH -A TG-PHY250035		# Project name
# #SBATCH -A TG-PHY230118		# Project name

# simulation parameters
t_end=2.5001
uemax=0.5

# flow variables
mul=5.00e-3
mug=2.50e-2
ul=0.5
ug=53.33 # so gas outlet velocity = 160 (Ao/Ai = 11.25)
theta0=70

export OMP_NUM_THREADS=48
./blur-ibm-openmpi $t_end $mul $mug $ul $ug $theta0 $uemax 1> out 2> log
#ibrun ./blur-ibm-openmpi $maxlevel $minlevel $tend $rhor $rhol $sigma $mul $mug $ul $ug $theta0 $uemax $tout $femax $semax $l0 1> out 2> log

