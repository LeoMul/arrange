# arrange
Arranges OMEGAXXXX into OMEGAZ. The motivation for the refactor is the old code stopping at say OMEGA0009 if OMEGA0010 does not exists - regardless if OMEGA0011 does. This was presumably for historical reasons. With large numbers of subworlds on current runs of the code I decided to refactor the code to dump all the data it can find as opposed to only sequentially labelled omega files.

While I was at it I updated it to the f90 standard. 

gfortran arrange.f90 -o arrange.x -O3

./arrange.x
