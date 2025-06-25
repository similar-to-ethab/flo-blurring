
# mpi
qcc -autolink -O2 -disable-dimensions -source -fopenmp blur-ibm.c

# hyperion
#mpicc -Wall -D_XOPEN_SOURCE=700 -std=c99 -O2 _blur-ibm.c -o blur3D -L$HOME/gl -lglutils -lfb_tiny -lm 

# local
gcc -Wall -D_XOPEN_SOURCE=700 -std=c99 -O2 _blur-ibm.c -o blur-ibm-openmp -L$BASILISK/gl -lglutils -lfb_tiny -lm 
