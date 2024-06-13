#!/bin/bash

# simulation_names=("sim3")

simulation_names=("sim4" "sim5" "sim6" "sim7") #"sim8" "sim9" "sim10" "sim11")
# simulation_names=("sim12" "sim13" "sim14" "sim15" "sim16" "sim17" "sim18" "sim19" "sim20")
deam=("no_deam" "deam")
coverage_values=("4")


# Loop through each combination of simulation name, deam value, and coverage value
for sim in "${simulation_names[@]}"; do
    for d in "${deam[@]}"; do
        for cov in "${coverage_values[@]}"; do
            # Construct the screen session name and log file path
            session_name="relate_${sim}_${d}_cov${cov}"
            log_file="${sim}/${d}/cov${cov}/many_kinships.log"

            # Create a new screen session with the specified name and run the script
            screen -S "$session_name" -L -Logfile "$log_file" -dm bash -c "source ./relate.sh $sim $d $cov"
        done
    done
done