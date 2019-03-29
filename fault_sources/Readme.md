## Code to reproduce Figures 23 to 28

Libraries required:

Install LBFGS library from https://github.com/chokkan/liblbfgs

Scripts to be run (in order) for a particular value of regularisation parameter c:

1. fault_sources.m : Matlab code to create a spatiotemporal source with several sources placed along the fault of the overthrust model

2. scons : Run SConstruct to compile the modelling and inversion codes

3. modeling_slurm.sh : SLURM script to perform modelling. Note: need to edit the path of the executable in this script

4. inversion_slurm.sh: SLURM script to perform inversion Note: need to edit the path of the executable in this script

5. blobslices.py : Python script to isolate blobs and plot them

6. plotdobsvel.py: Python script to plot data and velocity model

Main codes:

modeling_spatiotemporal.c : modelling code.  
Input->  win_overt.bin (Hard velocity model)  Output ->   dobs.bin (data)

lbfwi.c: inversion code*
Input -> dobs.bin (data) sm_win_overt.bin (Smooth velocity model)  Output -> inverted_src3d.bin

*Note: The sparsity promoting parameter c needs to be specified inside the code
