#!/bin/bash
#SBATCH --job-name=R_simulation_SLURM_attempt1       # Job name
#SBATCH --time=00:10:00  # Wall clock time limit
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=18
#SBATCH --gpus-per-node=1
#SBATCH --partition=gpu
#SBATCH --mail-type=ALL
#SBATCH --mail-user=l.s.b.yuan@uva.nl

# Load R module
module purge
module load 2023
module load R/4.3.2-gfbf-2023a

# echo "Start task $(date)" 
# echo -e "Directory: $(pwd) \n" 

# Create necessary directories in the scratch space
SCRATCH_DIR=/scratch-shared/$USER/${SLURM_JOB_ID}
mkdir -p $SCRATCH_DIR
cp /projects/0/prjs0958/SLURM_check/Simulation_downsampling_exit_time_right_deepening_iter1000.R $SCRATCH_DIR
cp /projects/0/prjs0958/SLURM_check/helper_functions_Rinn.R $SCRATCH_DIR

# Copy custom R libraries to the scratch space
CUSTOM_R_LIBS=/projects/0/prjs0958/R/library
SCRATCH_R_LIBS=$SCRATCH_DIR/R_libs
mkdir -p $SCRATCH_R_LIBS
cp -r $CUSTOM_R_LIBS/* $SCRATCH_R_LIBS

# Set R_LIBS_USER to point to the custom libraries in the scratch space
export R_LIBS_USER=$SCRATCH_R_LIBS

# Navigate to the scratch directory
cd $SCRATCH_DIR

# echo -e "Scratch Directory: $(pwd) \n" 

# Run the R script
Rscript Simulation_downsampling_exit_time_right_deepening_iter1000.R

# Copy output files back to the project directory
cp -r $SCRATCH_DIR/SLURM_downsample_figs /projects/0/prjs0958/SLURM_check/
cp -r $SCRATCH_DIR/SLURM_downsample_est /projects/0/prjs0958/SLURM_check/
cp $SCRATCH_DIR/makeForkCluster_log.txt /projects/0/prjs0958/SLURM_check/
cp $SCRATCH_DIR/execution_time_log.txt /projects/0/prjs0958/SLURM_check/

# Clean up
# rm -rf $SCRATCH_DIR

