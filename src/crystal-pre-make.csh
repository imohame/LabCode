#!/bin/csh
rm *.o *.mod
rm crystal2d
# dos2unix *
#csh -f /usr/local/apps/openmpi/ompi184_pgi151.csh
csh -f /usr/local/apps/openmpi/pgi151_ompi.csh
echo ---for pgi   ---     source /usr/local/apps/openmpi/pgi151_ompi.csh
echo ---for intel ---     source /usr/local/apps/openmpi/intel2013_ompi.csh
add pgi64_hydra
echo -----add pgi64_hydra
