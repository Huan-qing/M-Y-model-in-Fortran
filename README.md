The folder "model w gravitation" contains all the files needed to run the M-Y model (simulating particle transport considering the effect of gravitational sinking)

The folder "model w primary production" is a work in progress, as I tried to add the primary production to the M-Y model in addition to gravitational sinking (simulating particle transport considering the effect of gravitational sinking and primary production)

files:
mymodelXXX: the fortran file where the model code was saved 
obl: define the mixed layer depth as well as the types of parameters used in the model 
runpomXXX: the compiling file used to call the model in the Cygwin Terminal

The model was run in Cygwin 64 Terminal (https://www.cygwin.com/index.html)

Step1: get into the directory with runpom and all other files

$ cd ./directory name/

Step2: run the model

$ ./runpom
