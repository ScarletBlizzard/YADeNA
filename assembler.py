from argparse import ArgumentParser
from collections import defaultdict
import itertools

from Bio import SeqIO

import util


def create_parser():
    """Returns parser for command line arguments."""
    parser = ArgumentParser()
    # Required positional arguments
    parser.add_argument('reads_file',
                        help='File containing reads in FASTQ format')
    parser.add_argument('contigs_file',
                        help='File containing left and right contigs in FASTA format (left goes first)')
    parser.add_argument('alignments_file',
                        help='File containing read alignments')
    # Optional arguments
    parser.add_argument('-r', '--read_len', type=int, default=100,
                        help='Length of a read')
    parser.add_argument('-o', '--min_overlap_len', type=int,
                        help='Minimum length of reads overlap')
    parser.add_argument('-a', '--alignments_type', choices='art',
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
    """
    Receives reads, contigs and minimum overlap length.
    Returns overlap graph created in a naive way using hash table.
    """
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


class WrongPathError(Exception):
    """Raised when last read in path didn't match the read we want."""
    pass

def traverse(graph, reads, read, last, visited=set()):
    """
    Recursive function that traverses the overlap graph, thus assembling the
    target sequence. Uses greedy approach.

    Receives graph, current vertex (read), dict of reads, read that must be
    last in path and set of visited reads.
    Returns target sequence.
    """
    visited.add(read.id)
    descendants = sorted(graph[read.id], key=lambda x: x[1])
    while len(descendants) > 0:
        read_id, o_len = descendants.pop()
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


def main():
    """Returns target sequence."""
    parser = create_parser()
    args = parser.parse_args()

    reads = SeqIO.to_dict(SeqIO.parse(args.reads_file, 'fastq'))
    contigs = SeqIO.to_dict(SeqIO.parse(args.contigs_file, 'fasta'))
    if len(contigs) != 2:
        raise ValueError('contigs_file must contain only left and right contigs in FASTA format (left goes first)')
    if ids := set(reads) & set(contigs): # check for duplicate ids
        raise ValueError(f'Some reads and contigs have the same ids: {*ids,}')

    if args.alignments_type == 'art':
        orientation = util.parse_art_orientation(args.alignments_file)
        orient(reads, orientation)

    graph = create_overlap_graph(reads, contigs, args.min_overlap_len or args.read_len//15)

    (l_cont_id, l_cont), (r_cont_id, r_cont), *_ = contigs.items()
    reads[l_cont_id], reads[r_cont_id] = l_cont, r_cont # treat contigs like reads
    return traverse(graph, reads, l_cont, r_cont)


if __name__ == '__main__':
    print(main())
