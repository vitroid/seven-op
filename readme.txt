*** Compile ***

Eigen package is required to compile this program.
http://eigen.tuxfamily.org/index.php?title=Main_Page

In my case, eigen is installed via Mac HomeBrew package manager
% brew install eigen

and the program is compiled with specifying the location of header files explicitly.
g++  -O6 -I /usr/local/Cellar/eigen/3.1.1/include/eigen3 seven-op-gro.cpp -o seven-op-gro



*** usage ***

./seven-op-gro < coords.gro


*** benchmark ***

On my MBA, i runs more than 120x faster than the python version!
