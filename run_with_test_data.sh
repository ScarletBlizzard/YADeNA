#!/bin/bash

# This script runs assembler with data from dir containing test data and reports errors

max_tests_count=${1:-10} # How many times to run assembler with different data
read_len=${2:-100}
ins_len_mean=${3:-350}
test_data_dir=test_data
reads_file=${test_data_dir}/.reads.fq # temp file for storing reads
contigs_file=${test_data_dir}/.contigs.fa # temp file for storing contigs

i=1
for target_dir in ${test_data_dir}/out_dir/*; do
  # Getting target sequence
  target_file=${target_dir}/target.fq
  target_seq=$(cat ${target_file} | head -n 2 | tail -n 1)

  # Making reads file
  cat ${target_dir}/Sim_${read_len}_${ins_len_mean}_1.fq > $reads_file
  cat ${target_dir}/Sim_${read_len}_${ins_len_mean}_2.fq >> $reads_file

  # Making contigs file
  echo ">left_contig" > $contigs_file
  echo ${target_seq:0:read_len} >> $contigs_file # append left contig
  right_contig_pos="$((${#target_seq}-$read_len))"
  echo ">right_contig" >> $contigs_file
  echo ${target_seq:right_contig_pos:read_len} >> $contigs_file # append right contig

  # Running the assembler
  out=$(python3 assembler.py ${reads_file} ${contigs_file} -r ${read_len})

  # Checking and printing output
  if ! [[ $out == $target_seq ]]; then
    echo Unexpected output for ${target_dir}
  fi
  echo "$out"
  
  ((i++))
  if (( $i > $max_tests_count )); then
    break
  fi
done

rm $reads_file $contigs_file # remove temp files