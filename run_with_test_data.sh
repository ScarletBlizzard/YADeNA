#!/bin/bash

# This script runs assembler with data from dir containing test data and reports errors

max_tests_count=${1:-10} # How many times to run assembler with different data
read_len=${2:-150}
test_data_dir=test_data
contig_file_prefix=${test_data_dir}/.contig # prefix for temp files for storing contigs

# Workaround for running assembler.py with Python 3 both on Linux and Windows
python3_cmd=python
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  python3_cmd=python3
fi

i=1
for target_dir in ${test_data_dir}/out_dir/*; do
  # Getting target sequence
  target_file=${target_dir}/target.fa
  target_seq=$(cat ${target_file} | tail -n 1)

  # Making contigs file
  echo ">left_contig" > ${contig_file_prefix}_1.fa
  echo ${target_seq:0:read_len} >> ${contig_file_prefix}_1.fa # append left contig
  right_contig_pos="$((${#target_seq}-$read_len))"
  echo ">right_contig" > ${contig_file_prefix}_2.fa
  echo ${target_seq:right_contig_pos:read_len} >> ${contig_file_prefix}_2.fa # append right contig

  # Running the assembler
  out=$(${python3_cmd} assembler.py -r1 ${target_dir}/dat1.fq -r2 ${target_dir}/dat2.fq \
                                    -c1 ${contig_file_prefix}_1.fa -c2 ${contig_file_prefix}_2.fa \
                                    -a1 ${target_dir}/dat1.aln -a2 ${target_dir}/dat2.aln -at art)

  # Checking and printing output
  if [[ $out != $target_seq ]]; then
    echo Unexpected output for ${target_dir}:
  else
    echo Output for ${target_dir}:
  fi
  echo "$out"
  echo
  
  ((i++))
  if (( $i > $max_tests_count )); then
    break
  fi
done

rm ${contig_file_prefix}_1.fa ${contig_file_prefix}_2.fa # remove temp files