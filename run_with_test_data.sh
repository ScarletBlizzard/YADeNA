#!/bin/bash

# This script runs assembler with data from dir containing test data and reports errors

read_len=${1:-90}
ins_len_mean=${2:-170}
test_data_dir=test_data

for target_dir in ${test_data_dir}/out_dir/*; do
  target_file=${target_dir}/target.fq
  target_seq=$(cat ${target_file} | head -n 2 | tail -n 1)
  left_contig=${target_seq:0:read_len}
  left_contig_file=${test_data_dir}/left_contig.txt # temp file for storing left contig
  echo $left_contig > $left_contig_file
  right_contig_pos="$((${#target_seq}-$read_len))"
  right_contig=${target_seq:right_contig_pos:read_len}
  right_contig_file=${test_data_dir}/right_contig.txt # temp file for storing right contig
  echo $right_contig > $right_contig_file
  left_reads_file=${target_dir}/Sim_${read_len}_${ins_len_mean}_1.fq
  right_reads_file=${target_dir}/Sim_${read_len}_${ins_len_mean}_2.fq
  out=$(python3 assembler.py -lr ${left_reads_file} -rr ${right_reads_file} -lc ${left_contig_file} -rc ${right_contig_file})
  if ! [[ $out == $target_seq ]]; then
    echo Error on ${target_dir}
  fi
  echo $out
done
