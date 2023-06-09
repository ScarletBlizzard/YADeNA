"""
This module contains functions that are not quite related to the assembly
itself.
"""
from collections import Counter
from io import StringIO
import os
from pathlib import Path
import subprocess
import tempfile

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline as blastn
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from frozendict import frozendict
import pysam


Alignments = frozendict[str, pysam.AlignedSegment]


def compute_identity_of_unaligned(query_fname, subject_fname):
    output = blastn(
            query=query_fname, subject=subject_fname, outfmt=5)()[0]
    blast_record = NCBIXML.read(StringIO(output))
    for aln in blast_record.alignments:
        for hsp in aln.hsps:
            return hsp.identities / hsp.align_length
    return 0


class BAMUtil:
    def __init__(self, ref_fname: str, query_fnames: list[str], cleanup=True):
        with tempfile.NamedTemporaryFile(delete=False) as bam_file:
            try:
                subprocess.check_output(['bwa', 'index', ref_fname],
                                        stderr=subprocess.STDOUT)
                subprocess.check_output(['bwa', 'mem', ref_fname,
                                         '-o', bam_file.name, *query_fnames],
                                        stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                exit(e.output.decode('ascii'))
            self._bam_fname = bam_file.name
        pysam.sort('-o', self._bam_fname, self._bam_fname)
        pysam.index(self._bam_fname)

        self._ref_fname = ref_fname
        ref_record = next(SeqIO.parse(ref_fname, 'fasta'))
        self._ref_id = ref_record.id
        self._ref_seq = ref_record.seq
        self._cleanup = cleanup

    @classmethod
    def from_ref_str(cls, ref_seq: str, ref_id: str,
                     query_fnames: list[str], cleanup=True):
        if not ref_seq:
            raise ValueError('ref_seq must be a str of non-zero length')
        ref_fname = ''
        with tempfile.NamedTemporaryFile(delete=False, mode='w') as ref_file:
            SeqIO.write([SeqRecord(seq=Seq(ref_seq), id=ref_id)],
                        ref_file, 'fasta')
            ref_fname = ref_file.name
        return cls(ref_fname, query_fnames, cleanup)

    def save_sam(self):
        subprocess.run(['samtools', 'view', '-h',
                        '-o', 'output.sam', self._bam_fname])

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        Path(self._bam_fname).unlink(missing_ok=True)
        # Removing BAM index file and internal BWA files
        Path(self._bam_fname + '.bai').unlink(missing_ok=True)
        if self._cleanup:
            for ext in ('amb', 'ann', 'bwt', 'fai', 'pac', 'sa'):
                Path(f'{self._ref_fname}.{ext}').unlink(missing_ok=True)

    def consensus(self, start=0, end=None):
        if end is None:
            end = len(self._ref_seq)-1
        consensus_seq = ['']*len(self._ref_seq)

        with pysam.AlignmentFile(self._bam_fname, 'rb') as bam_file:
            nucleotide_counter = Counter()
            for pileup_column in bam_file.pileup(self._ref_id, start, end+1,
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
                    consensus_seq[pileup_column.reference_pos] = (max(
                            nucleotide_counter, key=nucleotide_counter.get))

        for i in range(len(consensus_seq)):
            if not consensus_seq[i]:
                consensus_seq[i] = self._ref_seq[i]

        return ''.join(consensus_seq)

    def fetch_alignments(self, start=0, end=None) -> Alignments:
        if end is None:
            end = len(self._ref_seq)

        alignments = {}
        with pysam.AlignmentFile(self._bam_fname, 'rb') as bam_file:
            for pileup_column in bam_file.pileup(self._ref_id, start, end+1):
                for pileup_read in pileup_column.pileups:
                    aln = pileup_read.alignment
                    read_id = aln.query_name
                    if read_id not in alignments:
                        alignments[read_id] = aln

        return frozendict(alignments)
