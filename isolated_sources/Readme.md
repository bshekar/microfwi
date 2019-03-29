## Code to reproduce Figures 17 and 18

Libraries required:

Install LBFGS library from https://github.com/chokkan/liblbfgs

Scripts to be run (in order) for a particular value of regularisation parameter c:

1. sources.py : Python script to generate four source wavelets

2. scons : Run SConstruct to compile the modelling and inversion codes

3. modeling_slurm.sh : SLURM script to run modelling code. Note: Edit the path in this script

4. inversion_slurm.sh: SLURM script to run inversion code. Note: Edit the path in this script

5. blobslices.py : Python script to isolate blobs and plot them. Outputs the isolated sources and the wavelets extracted from the centre of the blobs.

6. plotdobsvel.py: Python script to plot the data and velocity model


Main codes:

modeling.c : modelling code.  
Input->  win_overt.bin (Hard velocity model)  Output ->   dobs.bin (data)

lbfwi.c: inversion code*

Input -> dobs.bin (data) win_overt.bin (velocity model)  
source1.bin source2.bin source3.bin source4.bin (Four source wavelets)

Output -> inverted_src3d.bin

*Note: The sparsity promoting parameter c needs to be specified inside the code
