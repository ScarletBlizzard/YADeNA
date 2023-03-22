#!/bin/bash

# This script takes reference sequence and generates simulated data based on this sequence for testing the assembler

cd "$(dirname "$0")" # change to dir of this script

ref_file=ref.fa # FASTA/FASTQ file containing reference sequence
out_dir=out_dir
read_len=${1:-100} # length of a read
ins_len_mean=${2:-350}

PIRS_DIR=${PIRS_DIR:-"/usr/local/share/pirs"} # set path of dir containing pIRS executable if not set already
alias pirs="${PIRS_DIR}/pirs"
profiles="-B ${PIRS_DIR}/Profiles/Base-Calling_Profiles/humNew.PE100.matrix.gz -I ${PIRS_DIR}/Profiles/InDel_Profiles/phixv2.InDel.matrix -G ${PIRS_DIR}/Profiles/GC-depth_Profiles/humNew.gcdep_100.dat"

ref_id=$(cat ${ref_file} | head -n 1) # id of reference sequence
ref_seq=$(cat ${ref_file} | head -n 2 | tail -n 1) # reference sequence
ref_seq_len=${#ref_seq}
log_file=log.txt
err_file=err.txt
# gap_len is length of gap in sequence which is going to be filled by assembler
rm $log_file $err_file
for gap_len in `seq 200 100 1000`; do
  target_seq_len=$(($read_len*2 + $gap_len)) # length of resulting (target) sequence (two contigs and filled gap between them)
  max_pos=$(($ref_seq_len-$target_seq_len)) # maximum starting position of target sequence
  pos="$(shuf -i 0-"$max_pos" -n 1)" # random starting position of target sequence
  target_seq=${ref_seq:pos:target_seq_len}
  target_dir=${out_dir}/sim_${gap_len}
  mkdir -p $target_dir
  target_file=${target_dir}/target.fq # file containing target sequence
  echo ${ref_id}_${gap_len} > $target_file
  echo $target_seq >> $target_file
  pirs simulate ${target_file} -l ${read_len} -x 5 -m ${ins_len_mean} ${profiles} -o ${target_dir} >> ${log_file} 2>>${err_file}
done
