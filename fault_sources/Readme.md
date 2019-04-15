## Code to reproduce Figures 23 to 28

Libraries required:

Install LBFGS library from http://www.chokkan.org/software/liblbfgs/ (I recommend downloading from the link labelled "source code." The library labelled from their github repo did not compile in my laptop.)

Scripts to be run (in order) for a particular value of regularisation parameter c:

1. fault_sources.m : Matlab code to create a spatiotemporal source with several sources placed along the fault of the overthrust model

2. scons : Run SConstruct to compile the modelling and inversion codes

3. ./modeling_spatiotemporal.exe: Run the modeling code with OpenMP and all available threads.

3a. modeling_slurm.sh : SLURM script to run modelling code. Note: Use this ONLY if you have SLURM scheduler installed
Sample execution time: 8 seconds on  Intel Xeon Processor E5-2680 v3 24 Cores 2.5GHz 30MB 2133MHz

4. ./lbfwi.exe: Run the inversion code for 4 iterations with OpenMP and all available threads.

4a. inversion_slurm.sh : SLURM script to run modelling code. Note: Use this ONLY if you have SLURM scheduler installed
Sample execution time for 4 iterations: 290 seconds on Intel Xeon Processor E5-2680 v3 24 Cores 2.5GHz 30MB 2133MHz

5. blobslices.py : Python script to isolate blobs and plot them

6. plotdobsvel.py: Python script to plot data and velocity model

Main codes:

modeling_spatiotemporal.c : modelling code.  
Input->  win_overt.bin (Hard velocity model)  Output ->   dobs.bin (data)

lbfwi.c: inversion code*
Input -> dobs.bin (data) sm_win_overt.bin (Smooth velocity model)  Output -> inverted_src3d.bin

*Note: The sparsity promoting parameter c and number of iterations need to be specified inside the code
