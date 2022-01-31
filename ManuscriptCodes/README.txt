These folders contains all of the files necessary to run this code. I have created a make file that gathers all of the necessary files to run the main scripts in each folder (DDC_2D_mpi.f90 and DDC_3D_mpi.f90). 

Each makefile creates an executable that can be submitted for a job on a cluster with a slurm file or can be run locally on fewer numbers of cores. Both folders also contain tiling algorithms for users to input simulation parameters and receive suggested tilings based on a desired efficiency level.

The results from the tiling programs will need to be updated in the main scripts for the variables nnx, nny, and nnz (if using the 3-d program) before the job is submitted. 

Hopefully the files are annotated such that they are easily read, but please reach out to me at lschauer@mines.edu if you have any questions. 
