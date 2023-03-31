#!/bin/bash

# This script runs assembler with data from dir containing test data and reports errors

max_tests_cnt=${1:-10} # How many times to run assembler with different data
read_len=${2:-150}
test_data_dir=test_data
reads_file=${test_data_dir}/.reads.fq # temp file for storing reads
contigs_file=${test_data_dir}/.contigs.fa # temp file for storing contigs
alignments_file=${test_data_dir}/.alignments.aln # temp file for storing contigs

log=test_runs.log
rm $log

# Workaround for running assembler.py with Python 3 both on Linux and Windows
python3_cmd=python
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  python3_cmd=python3
fi

ran_tests_cnt=0
successful_cnt=$max_tests_cnt
for target_dir in ${test_data_dir}/out_dir/*; do
  # Getting target sequence
  target_file=${target_dir}/target.fa
  target_seq=$(cat ${target_file} | head -n 2 | tail -n 1)

  # Making reads file
  cat ${target_dir}/dat1.fq > $reads_file
  cat ${target_dir}/dat2.fq >> $reads_file

  # Making contigs file
  echo ">left_contig" > $contigs_file
  echo ${target_seq:0:read_len} >> $contigs_file # append left contig
  right_contig_pos="$((${#target_seq}-$read_len))"
  echo ">right_contig" >> $contigs_file
  echo ${target_seq:right_contig_pos:read_len} >> $contigs_file # append right contig

  # Making alignments file
  grep -o '^>.*' ${target_dir}/dat1.aln > $alignments_file
  grep -o '^>.*' ${target_dir}/dat2.aln >> $alignments_file

  echo [${target_dir}] >> $log

  # Running the assembler
  out=$(${python3_cmd} assembler.py ${reads_file} ${contigs_file} ${alignments_file} -r ${read_len} -a art 2>>${log})

  if [[ $out != $target_seq ]]; then
    ((successful_cnt--))
    echo Unexpected result: >> $log
    echo "$out" >> $log
    echo Expected: >> $log
    echo $target_seq >> $log
  else
    echo OK >> $log
  fi
  echo >> $log
  
  ((ran_tests_cnt++))
  if (( $ran_tests_cnt > $max_tests_cnt )); then
    break
  fi
done

echo ${successful_cnt}/${max_tests_cnt} tests ran successfully

rm $reads_file $contigs_file $alignments_file # remove temp files