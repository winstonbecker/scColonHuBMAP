#################
#set a job name  
#SBATCH --job-name=j2i2
#################  
#a file for job output, you can check job progress
#SBATCH --output=j2i2.out.txt
#################
# a file for errors from the job
#SBATCH --error=j2i2.err
#################
#time you think you need; default is one hour
#in minutes in this case, hh:mm:ss
#SBATCH --time=1:00:00
#################
#quality of service; think of it as job priority
#SBATCH --partition=normal
#################
#number of nodes you are requesting
#SBATCH --mem=64G
#################
ml biology
ml R/4.1.2

# Define paths
i_var=$1
j_var=$2

srun Rscript ./2b_run_decontX.R $i_var $j_var