"""
This module contains functions that are not quite related to the assembly
itself.
"""
from collections import Counter
from io import StringIO
import os
import subprocess
import tempfile

from Bio.Blast.Applications import NcbiblastnCommandline as blastn
from Bio.Blast import NCBIXML
from frozendict import frozendict
import pysam


def compute_identity_of_unaligned(query_fname, subject_fname):
    output = blastn(
            query=query_fname, subject=subject_fname, outfmt=5)()[0]
    blast_record = NCBIXML.read(StringIO(output))
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            return hsp.identities / hsp.align_length
    return 0


def compute_identity(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError('Sequences must be aligned and have the same length')
    match_count = 0
    alignment_len = len(seq1)
    for i in range(alignment_len):
        if seq1[i] == seq2[i]:
            match_count += 1
    return match_count / alignment_len


def align_SAM(sam_fname, ref_fname, query_fnames):
    """
    Name of sorted SAM file to create, name of FASTA/FASTQ file
    containing reference sequence and names of FASTA/FASTQ files
    containing query sequence which to map onto the
    given sequence using minimap2 aligner.
    """
    with open(sam_fname, 'w', encoding='utf-8') as sam_file:
        subprocess.run(['minimap2', '--MD', '-a', ref_fname, *query_fnames],
                       stdout=sam_file,
                       stderr=subprocess.DEVNULL,
                       check=True)
    pysam.sort('-o', sam_fname, sam_fname)


class BAMUtil:
    def __init__(self, ref_fname, ref_id, query_fnames):
        """
        Name of sorted BAM file to create and index, name of FASTA/FASTQ file
        containing reference sequence and names of FASTA/FASTQ files
        containing query sequence which to map onto the
        given sequence using BWA.
        """
        with tempfile.NamedTemporaryFile(delete=False) as bam_file:
            subprocess.run(['bwa', 'index', ref_fname],
                           stderr=subprocess.DEVNULL,
                           check=True)
            subprocess.run(['bwa', 'mem', ref_fname, *query_fnames],
                           stdout=bam_file,
                           stderr=subprocess.DEVNULL,
                           check=True)
            self._bam_fname = bam_file.name
        pysam.sort('-o', self._bam_fname, self._bam_fname)
        pysam.index(self._bam_fname)
        self._ref_id = ref_id

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        os.remove(self._bam_fname)
        os.remove(self._bam_fname + '.bai')  # remove index file

    def save(self, bam_fname):
        pass

    @staticmethod
    def _parse_cigar(seq, cigartuples):
        prepared_seq = []
        idx = 0
        for op, length in cigartuples:
            if op == 0:  # M
                prepared_seq.append(seq[idx:idx+length])
                idx += length
            elif op == 2:  # D
                idx += length
            elif op in (1, 4):  # I or S
                prepared_seq.append('-'*length)
            else:
                print(op)
        return ''.join(prepared_seq).upper()

    def compute_avg_identity(self, start, stop):
        identity_sum = 0
        alignments_count = 0
        with pysam.AlignmentFile(self._bam_fname, 'rb') as bam_file:
            for alignment in bam_file.fetch(self._ref_id, start, stop+1):
                alignments_count += 1
                if alignment.cigartuples:
                    ref_part = self._parse_cigar(
                            alignment.get_reference_sequence(),
                            alignment.cigartuples)
                    identity_sum += compute_identity(
                            alignment.query_sequence, ref_part)

        if alignments_count > 0:
            return identity_sum / alignments_count
        else:
            return 0

    def consensus(self, start, stop):
        """
        Receives name of BAM file containing assembled sequence and reads
        mapped to it.
        Returns their consensus sequence.
        """
        with pysam.AlignmentFile(self._bam_fname, 'rb') as bam_file:
            consensus_seq = []
            nucleotide_counter = Counter()
            for pileup_column in bam_file.pileup(self._ref_id, start, stop+1,
                                                 truncate=True):
                nucleotide_counter.clear()

                for pileup_read in pileup_column.pileups:
                    try:
                        nucleotide = pileup_read.alignment.query_sequence[
                                pileup_read.query_position]
                        nucleotide_counter[nucleotide] += 1
                    except TypeError:
                        pass

                if nucleotide_counter:
                    consensus_seq.append(max(nucleotide_counter,
                                             key=nucleotide_counter.get))

        return ''.join(consensus_seq)

    def fetch_reads(self, start, stop):
        reads, read_positions = {}, {}
        with pysam.AlignmentFile(self._bam_fname, 'rb') as bam_file:
            for pileup_column in bam_file.pileup(self._ref_id, start, stop+1,
                                                 truncate=False):
                for pileup_read in pileup_column.pileups:
                    alignment = pileup_read.alignment
                    read_id = alignment.query_name
                    if read_id not in reads:
                        reads[read_id] = alignment.query_sequence
                        read_positions[read_id] = alignment.reference_start
        return frozendict(reads), read_positions
