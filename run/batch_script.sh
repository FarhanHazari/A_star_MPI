#!/bin/bash

# Init color to print messages
C_RESET='\033[0m'
C_PURPLE='\033[0;35m'

# infoln echos in purple color
function infoln() {
  echo -e "${C_PURPLE}${1}${C_RESET}"
}

seed=99
dim_start=256
type=maze
algo=1

current_dir=$(dirname "$(realpath -s "$0")")
project_dir=$(dirname "${current_dir}")
exec_file="${project_dir}/a_star"

#Clean and build executable
infoln "Building executable.."
module load mpich-3.2
cd $project_dir && make clean && make

outputs_dir="${current_dir}/outputs"
errors_dir="${current_dir}/errors"
scripts_dir="${current_dir}/scripts"

#Delete previous stuff
infoln "Removing previous outputs if any.."
rm -f ${outputs_dir}/* ${errors_dir}/* ${scripts_dir}/*

#Create directories
mkdir -p $outputs_dir
mkdir -p $errors_dir
mkdir -p $scripts_dir

#Loop to set the number of cores
for c in {0..6}
do
nb_cores=$((2**c))

#Loop to set dimensions
for d in {0..4}
do
dim=$((dim_start*2**$d))

name=n${nb_cores}-d${dim}
script_file=${name}.sh
script_path=${scripts_dir}/${script_file}

#Create the script file
infoln "Creating script for ${nb_cores} cores and a ${dim}x${dim} grid"
cat <<EOF > $script_path
#!/bin/bash

#PBS -N $name
#PBS -l select=${nb_cores}:ncpus=1:mem=2gb -l place=pack:excl
#PBS -l walltime=0:30:00
#PBS -q short_cpuQ

module load mpich-3.2
mpirun.actual -n $nb_cores $exec_file $seed $dim $dim $type $algo
EOF

#Make the file executable
chmod +x $script_path

infoln "Executing script 10 times"

for i in {1..10}
do

output_file=${outputs_dir}/${name}-${i}-o.txt
error_file=${errors_dir}/${name}-${i}-e.txt

#Run the qsub command
qsub -o $output_file -e $error_file $script_path 

#Check that the job has been submitted and we haven't reached the queue limit
while [ $? -ne 0 ]
do
    #If limit has been reached, sleep 30s and try again      
    echo "We have reached the queue limit.. Sleeping for 30s"
    sleep 30s
    qsub -o $output_file -e $error_file $script_path 
done

done

done

done