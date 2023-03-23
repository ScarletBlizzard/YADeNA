from argparse import ArgumentParser
from collections import defaultdict
import itertools

from Bio import SeqIO


def create_parser():
    """Returns parser for command line arguments."""
    parser = ArgumentParser()
    # Required positional arguments
    parser.add_argument('reads_file',
                        help='File containing reads in FASTQ format')
    parser.add_argument('contigs_file',
                        help='File containing left and right contigs in FASTA format (left goes first)')
    # Optional arguments
    parser.add_argument('-r', '--read_len', type=int, default=100,
                        help='Length of a read')
    parser.add_argument('-o', '--min_overlap_len', type=int,
                        help='Minimum length of reads overlap')
    return parser


def overlap_len(s1, s2, min_len):
    """
    Receives strings s1 and s2.
    Returns the length of the longest suffix of s1 that matches the prefix of s2.
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
    """Creates overlap graph in a naive way and returns it."""
    # TODO: Use smarter data structure and create graph in a smarter way
    graph = defaultdict(list)
    # Find overlaps between different reads, including between paired-end reads
    for r1, r2 in itertools.permutations(reads.values(), 2):
        if r1.seq != r2.seq:
            o_len = overlap_len(r1.seq, r2.seq, min_overlap_len)
            if o_len > 0:
                graph[r1.id].append((r2.id, o_len))
                
    c1, c2 = contigs.values()
    for r in reads.values():
        # Find overlaps between suffix of left contig and prefixes of all reads
        o_len = overlap_len(c1.seq, r.seq, min_overlap_len)
        if o_len > 0:
            graph[c1.id].append((r.id, o_len))
        # Find overlaps between prefix of right contig and suffixes of all reads
        o_len = overlap_len(r.seq, c2.seq, min_overlap_len)
        if o_len > 0:
            graph[r.id].append((c2.id, o_len))
    return graph


def main():
    parser = create_parser()
    args = parser.parse_args()
    # TODO: Raise exception if there are duplicate ids
    reads = SeqIO.to_dict(SeqIO.parse(args.reads_file, "fastq"))
    contigs = SeqIO.to_dict(SeqIO.parse(args.contigs_file, "fasta"))
    graph = create_overlap_graph(reads, contigs, args.min_overlap_len or args.read_len//2)
    print(graph)
    # TODO: Traverse graph to get the sequence


if __name__ == '__main__':
    main()
