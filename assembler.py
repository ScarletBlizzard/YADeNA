from argparse import ArgumentParser
import itertools
import os
import subprocess

from Bio import SeqIO

from overlap_graph import OverlapGraph
import util


def create_parser():
    """Returns parser for command line arguments."""
    parser = ArgumentParser()
    # Required arguments
    required = parser.add_argument_group('required named arguments')
    required.add_argument('-r', '--read_len', required=True, type=int,
                          help='Length of a read')
    required.add_argument('-r1', '--reads1', required=True,
                          help='File containing reads in FASTQ format')
    required.add_argument('-r2', '--reads2', required=True,
                          help='File containing reads in FASTQ format')
    required.add_argument('-c1', '--contig1', required=True,
                          help='File containing left contig in FASTA format')
    required.add_argument('-c2', '--contig2', required=True,
                          help='File containing right contig in FASTA format')
    required.add_argument('-a1', '--alignments1', required=True,
                          help='File containing read alignments')
    required.add_argument('-a2', '--alignments2', required=True,
                          help='File containing read alignments')
    # Optional arguments
    parser.add_argument('-o', '--output',
                        help='File containing fully assembled sequence')

    parser.add_argument('-d', '--read_len_divisor', type=int, default=3,
                        help=('How many times the minimum overlap length is '
                              'smaller than the read length '
                              '(used when min_overlap_len not given)'))
    parser.add_argument('-ol', '--min_overlap_len', type=int,
                        help='Minimum length of reads overlap')
    parser.add_argument('-at', '--alignments_type', choices=['art'],
                        help='Name of program that generated alignments files')
    return parser


def overlap_len(s1, s2, min_len):
    """
    Receives strings s1 and s2.
    Returns the length of the longest suffix of s1 that matches
    the prefix of s2.
    """
    start = 0
    while True:
        start = s1.find(s2[:min_len], start)
        if start == -1:
            return 0
        if s2.startswith(s1[start:]):
            return len(s1)-start
        start += 1


def create_overlap_graph(reads, contigs, min_overlap_len):
    """
    Receives reads, contigs and minimum overlap length.
    Returns overlap graph created in a naive way using hash table.
    """
    # TODO: Use smarter data structure and create graph in a smarter way
    graph = OverlapGraph()
    # Find overlaps between different reads, including between paired-end reads
    for r1, r2 in itertools.permutations(reads.values(), 2):
        if r1.seq != r2.seq:
            o_len = overlap_len(r1.seq, r2.seq, min_overlap_len)
            if o_len > 0:
                graph.add_child(r1.id, r2.id, o_len)

    c1, c2 = contigs.values()
    for r in reads.values():
        # Find overlaps between suffix of left contig and prefixes of reads
        o_len = overlap_len(c1.seq, r.seq, min_overlap_len)
        if o_len > 0:
            graph.add_child(c1.id, r.id, o_len)
        # Find overlaps between prefix of right contig and suffixes of reads
        o_len = overlap_len(r.seq, c2.seq, min_overlap_len)
        if o_len > 0:
            graph.add_child(r.id, c2.id, o_len)
    return graph


class WrongPathError(Exception):
    """Raised when last read in path didn't match the read we want."""


def traverse(graph, reads, read, last, visited=None):
    """
    Recursive function that traverses the overlap graph, thus assembling the
    target sequence. Uses greedy approach.

    Receives graph, current vertex (read), dict of reads, read that must be
    last in path and set of visited reads.
    Returns target sequence.
    """
    if visited is None:
        visited = set()
    visited.add(read.id)
    while child := graph.get_next_child(read.id):
        read_id, o_len = child

        if read_id in visited:
            continue
        try:
            res = traverse(graph, reads, reads[read_id], last, visited)
            return str(read.seq) + res[o_len:]
        except WrongPathError:
            continue
    if read.id != last.id:
        raise WrongPathError
    return str(last.seq)


def orient(reads, orientation):
    """
    Receives dict of reads and dict that maps their ids to orientations
    which are either '+' or '-'. '-' means that the read is reverse complement
    of original sequence.
    Modifies reads dict according to orientation.
    """
    for r_id, r in reads.items():
        if orientation[r_id] == '-':
            r.seq = r.seq.reverse_complement()


def assemble(
        reads, contigs, orientation, read_len,
        read_len_divisor, min_overlap_len=None):
    """Returns assembled pre-consensus sequence based on args."""
    if len(contigs) != 2:
        raise ValueError(('contigs_file must contain only left and right'
                          'contigs in FASTA format (left goes first)'))
    if ids := set(reads) & set(contigs):  # check for duplicate ids
        raise ValueError(f'Some reads and contigs have the same ids: {*ids,}')

    if orientation:
        orient(reads, orientation)

    (l_con_id, l_con), (r_con_id, r_con), *_ = contigs.items()
    reads[l_con_id], reads[r_con_id] = l_con, r_con  # treat contigs like reads

    min_overlap_len = min_overlap_len or read_len//read_len_divisor
    while True:
        graph = create_overlap_graph(reads, contigs, min_overlap_len)
        try:
            seq = traverse(graph, reads, l_con, r_con)
            return seq
        except WrongPathError:
            min_overlap_len -= 1


def run(args):
    """Returns fully assembled sequence based on args."""
    reads = {**SeqIO.to_dict(SeqIO.parse(args['reads1'], 'fastq')),
             **SeqIO.to_dict(SeqIO.parse(args['reads2'], 'fastq'))}
    contigs = {**SeqIO.to_dict(SeqIO.parse(args['contig1'], 'fasta')),
               **SeqIO.to_dict(SeqIO.parse(args['contig2'], 'fasta'))}

    orientation = None
    if args['alignments_type'] == 'art':
        orientation = util.parse_art_orientation(
                (args['alignments1'], args['alignments2']))

    seq = assemble(
            reads, contigs, orientation, args['read_len'],
            args['read_len_divisor'], min_overlap_len=args['min_overlap_len'])

    pre_consensus_fname = 'pre_consensus.fa'
    with open(pre_consensus_fname, 'w', encoding='utf-8') as file:
        file.write('>pre_consensus\n')
        file.write(seq)

    temp_sam_fname = 'temp.sam'
    with open(temp_sam_fname, 'w', encoding='utf-8') as sam_file:
        subprocess.run(['minimap2', '-a', pre_consensus_fname, args['reads1'],
                        args['reads2'], args['contig1'], args['contig2']
                        ],
                       stdout=sam_file, stderr=subprocess.DEVNULL, check=True)
        os.remove(pre_consensus_fname)

    seq = util.consensus(temp_sam_fname)
    os.remove(temp_sam_fname)
    return seq


def main():
    parser = create_parser()
    args = vars(parser.parse_args())
    seq = run(args)
    print(seq)
    if args['output']:
        with open(args['output'], 'w', encoding='utf-8') as file:
            file.write('>output\n')
            file.write(seq)


if __name__ == '__main__':
    main()
