#!/bin/bash

# # # Set the simulation names, deam values, and coverage values
# # simulation_names=("sim4" "sim5" "sim6" "sim7" "sim8" "sim9" "sim10" "sim11")
# simulation_names=("sim12" "sim13" "sim14" "sim15" "sim16" "sim17" "sim18" "sim19" "sim20")

# # simulation_names=("sim1" "sim2")
# # # deam=("no_deam" "deam")
# # # coverage_values=(1 2 4 8)

# # simulation_names=("sim12" "sim13" "sim14" "sim15" "sim16" "sim17")

# # simulation_names=("sim12")

# deam=("no_deam" "deam")
# coverage_values=(1 2 4 10)


# # Loop through each combination of simulation name, deam value, and coverage value
# for sim in "${simulation_names[@]}"; do
#     for d in "${deam[@]}"; do
#         for cov in "${coverage_values[@]}"; do
#             # Construct the screen session name and log file path
#             session_name="${sim}_${d}_cov${cov}"
#             log_file="${sim}/${d}/cov${cov}/many_kinships.log"

#             # Create a new screen session with the specified name and run the script
#             screen -S "$session_name" -L -Logfile "$log_file" -dm bash -c "source ./fullSim2.sh $sim $d $cov"
#         done
#     done
# done




# # # Set the simulation names, deam values, and coverage values
# simulation_names=("sim12" "sim13" "sim14")
# # deam=("no_deam" "deam")
# # coverage_values=(0.5)

# simulation_names=("sim20")

# # # dandy3
# simulation_names=(  "sim3" "sim9" "sim15" "sim16")
# # # dandy2
# simulation_names=("sim10" "sim11" "sim12" "sim13" "sim14" )
# # # #dandy1
# simulation_names=("sim15" "sim16" "sim17" "sim18" "sim19" "sim20")
# simulation_names=("sim1" "sim2" "sim3" "sim4" "sim5" "sim6")

# simulation_names=( "sim3" "sim4" "sim5" "sim6" "sim7" "sim8" "sim9" "sim10" "sim11" "sim12" "sim13" "sim14")
# simulation_names=("sim12" "sim13" "sim14" "sim15" "sim16" "sim17" )
# simulation_names=("sim4" "sim5" "sim6" "sim7" "sim8" "sim9")
# simulation_names=("sim1" "sim2" "sim3" "sim4" "sim5" "sim6" "sim7" "sim8" "sim9" "sim10" "sim11" "sim12" "sim13" "sim14")
simulation_names=( "sim15" "sim16" "sim17" "sim18" "sim19" "sim20")
deam=("no_deam" "deam")
coverage_values=(0.1 0.5 1 2 4 8 10 20)
# coverage_values=(10)


# Loop through each combination of simulation name, deam value, and coverage value
for sim in "${simulation_names[@]}"; do
    for d in "${deam[@]}"; do
        for cov in "${coverage_values[@]}"; do
            # Construct the screen session name and log file path
            session_name="correctKin_${sim}_${d}_cov${cov}"
            log_file="${sim}/${d}/cov${cov}/many_kinships.log"

            # Create a new screen session with the specified name and run the script
            screen -S "$session_name" -L -Logfile "$log_file" -dm bash -c "source ./read.sh $sim $d $cov "
        done
    done
done

            # screen -S "$session_name" -L -Logfile "$log_file" -dm bash -c "source ./fullSim2.sh $sim/$d/cov$cov allbams_final_KING_bf KING $PWD "


# ll=("founder0" "founder1" "founder10" "founder11" "founder12" "founder13" "founder14" "founder15" "founder16" "founder17" "founder17" "founder18" "founder19" "founder2" "founder20" "founder21" "founder22" "founder23" "founder24" "founder25" "founder26" "founder27" "founder28" "founder29" "founder3" "founder4" "founder5" "founder6" "founder7" "founder8" "founder9")         
# # for f in $(ls $fq_path/*.gz|cut -f1 -d.|sort -u) 

# for f in "founder0" "founder1" "founder10" "founder11" "founder12" "founder13" "founder14" "founder15" "founder16" "founder17" "founder17" "founder18" "founder19" "founder2" "founder20" "founder21" "founder22" "founder23" "founder24" "founder25" "founder26" "founder27" "founder28" "founder29" "founder3" "founder4" "founder5" "founder6" "founder7" "founder8" "founder9"
# do
#     bwa mem fasta/new_reference_genome.fa <(cat fq/$f.h1.final.trimm.fq.gz fq/$f.h2.final.trimm.fq.gz) -t 16 |samtools sort -@4 -m4G - > $bam_path/$f.bam
# done


# "founder0" "founder1" "founder10" "founder11" "founder12" "founder13" "founder14" "founder15" "founder16" "founder17" "founder18" "founder19" "founder2" "founder20" "founder21" "founder22" "founder23" "founder24" "founder25" "founder26" "founder27" "founder28" "founder29" "founder3" "founder4" "founder5" "founder6" "founder7" "founder8" "founder9"

# "founder13" "founder14" "founder15" "founder16" "founder17"  "founder18" "founder19" "founder2" "founder20" "founder21" "founder22" "founder23" "founder24" "founder25" "founder26" "founder27" "founder28" "founder29" "founder3" "founder4" "founder5" "founder6" "founder7" "founder8" "founder9"


# for f in  "founder21" "founder22" "founder24" "founder25" "founder26" "founder27" ;do bwa mem /projects/korneliussen/people/vlr112/final_simulations/fasta/new_reference_genome.fa <(cat fq/$f.h1.final.trimm.fq.gz fq/$f.h2.final.trimm.fq.gz) -t 16 |samtools sort -@4 -m4G - > bam/$f.bam ; done



# lcmlkin -i /projects/korneliussen/people/vlr112/final_simulations/sim10/deam/cov2/results/lcmlkin/allbams_mp_update.vcf  -o lcmlkin_results_update -g all


## kill all my screens in server
# for session in $(screen -ls | grep "sim"| awk '{print $1}' | cut -d. -f 2); do screen -S "${session}" -X quit; done 

# for session in $(screen -ls | grep "pts"| awk '{print $1}' | cut -d. -f 2); do echo $session ; done 


# # Construct the screen session name and log file path
# session_name="lcmlkin_faroe"
# log_file="faroe/many_kinships2.log"

# # Create a new screen session with the specified name and run the script
# screen -S "$session_name" -L -Logfile "$log_file" -dm bash -c "source ./correctKin_realData.sh "






