#!/bin/bash

# This script takes reference sequence and generates simulated data based on this sequence for testing the assembler

ref_file=ref.fa # FASTA/FASTQ file containing reference sequence
out_dir=out_dir
read_len=90 # length of a read

PIRS_DIR=${PIRS_DIR:-"/usr/local/share/pirs"} # set path of dir containing pIRS executable if not set already
alias pirs="${PIRS_DIR}/pirs"
profiles="-B ${PIRS_DIR}/Profiles/Base-Calling_Profiles/humNew.PE100.matrix.gz -I ${PIRS_DIR}/Profiles/InDel_Profiles/phixv2.InDel.matrix -G ${PIRS_DIR}/Profiles/GC-depth_Profiles/humNew.gcdep_100.dat"

ref_id=$(cat ${ref_file} | head -n 1) # id of reference sequence
ref_seq=$(cat ${ref_file} | head -n 2 | tail -n 1) # reference sequence
ref_seq_len=${#ref_seq}
# gap_len is length of gap in sequence which is going to be filled by assembler
for gap_len in `seq 200 100 1000`; do
  target_seq_len=$(($read_len*2 + $gap_len)) # length of resulting (target) sequence (two contigs and filled gap between them)
  max_pos=$(($ref_seq_len-$target_seq_len)) # maximal starting position of target sequence
  pos="$(shuf -i 0-"$max_pos" -n 1)" # random starting position of target sequence
  target_seq=${ref_seq:pos:target_seq_len}
  target_file=${out_dir}/target_${gap_len}.txt # file containing target sequence
  echo $ref_id $gap_len > $target_file
  echo $target_seq >> $target_file
  pirs simulate ${target_file} -l ${read_len} -x 5 -m 170 -v 10 -z ${profiles} -o ${out_dir}/sim_${gap_len} >> log.txt 2>>err.txt
done
