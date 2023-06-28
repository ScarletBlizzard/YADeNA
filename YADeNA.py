from argparse import ArgumentParser

from Bio import SeqIO

from assembler import Assembler
from util import BAMUtil, Alignments


def create_parser():
    """Returns parser for command line arguments."""
    parser = ArgumentParser()
    # Required arguments
    required = parser.add_argument_group('required named arguments')
    required.add_argument('-r', '--reference', required=True,
                          help='Reference sequence')
    required.add_argument('-l', '--read_len', required=True, type=int,
                          help='Length of a read')
    required.add_argument('-r1', '--reads1', required=True,
                          help='File containing reads in FASTQ format')
    required.add_argument('-r2', '--reads2', required=True,
                          help='File containing reads in FASTQ format')
    # TODO: Constrain minid to the range [0, 100]
    required.add_argument('-i', '--minid', required=True, type=float,
                          help=('Minimum average percent-identity of a '
                                'genome region to consider it a gap'))
    # Optional arguments
    parser.add_argument('-o', '--output',
                        help='File containing fully assembled sequence')

    parser.add_argument('-d', '--read_len_div', type=int, default=3,
                        help=('How many times the minimum overlap length is '
                              'smaller than the read length '
                              '(used when min_overlap_len not given)'))
    parser.add_argument('-ol', '--min_overlap_len', type=int,
                        help='Minimum length of reads overlap')
    parser.add_argument('-gv', '--visualize',
                        dest='visualize', action='store_true',
                        help='Visualize overlap graph using Graphviz')
    parser.set_defaults(visualize=False)
    return parser


def compute_identity(read, seq_part, cigartuples):
    match_count, length = 0, 0
    i, j, d = 0, 0, 0
    for op, ln in cigartuples:
        length += ln
        if op == 0:  # M
            for _ in range(ln):
                if read[i+d] == seq_part[j]:
                    match_count += 1
                i += 1
                j += 1
        else:
            i += ln
            if op == 2:  # D
                d += ln
    return match_count / length


def compute_avg_identity(bam_util, seq, start, end):
    identity_sum = 0
    alignments = bam_util.fetch_alignments(start, end)
    alignments_count = len(alignments)
    for aln in alignments.values():
        if aln.cigartuples:  # TODO: Consider unmapped reads too
            seq_part = seq[aln.reference_start:aln.reference_end]
            identity_sum += compute_identity(
                    aln.query_sequence, seq_part, aln.cigartuples)
    if alignments_count == 0:
        return 0
    return identity_sum / alignments_count


def find_gaps(bam_util, seq, read_len, min_avg_identity):
    gaps = []
    window_len = 2 * read_len
    start, step = 0, read_len
    while start + step < len(seq):
        end = min(start+window_len-1, len(seq)-1)
        avg_identity = compute_avg_identity(bam_util, seq, start, end)
        if avg_identity < min_avg_identity:
            if not gaps or start-gaps[-1][1] > step+1:
                gaps.append((start, end))
            else:
                prev_start, _ = gaps.pop()
                gaps.append((prev_start, end))
        start += step
    return gaps


def get_contigs(gaps, consensus_seq):
    if gaps:
        contigs = []
        if gaps[0][0] != 0:  # starting with a contig
            contigs.append(consensus_seq[0:gaps[0][0]])
        for i in range(1, len(gaps)):
            contigs.append(consensus_seq[gaps[i-1][1]+1:gaps[i][0]])
        if gaps[-1][1] != len(consensus_seq)-1:  # ending with a contig
            contigs.append(consensus_seq[gaps[-1][1]+1:len(consensus_seq)])
        return contigs
    return [consensus_seq]


def get_reads_with_positions(alignments: Alignments):
    reads = {read_id: aln.query_sequence
             for read_id, aln in alignments.items()}
    read_positions = {read_id: aln.reference_start for
                      read_id, aln in alignments.items()}
    return reads, read_positions


def main(args):
    ref_record = next(SeqIO.parse(args['reference'], 'fasta'))
    ref_id, ref_len = ref_record.name, len(ref_record.seq)
    read_len = args['read_len']
    with BAMUtil(args['reference'], [args['reads1'],
                 args['reads2']]) as bam_util:
        consensus_seq = bam_util.consensus()
        gaps = find_gaps(
                bam_util, consensus_seq, read_len, args['minid']/100)

        contigs = get_contigs(gaps, consensus_seq)
        if not gaps:
            return contigs[0]
        if not contigs:
            raise ValueError(
                    ('Cannot assemble sequence with given '
                     f'minimum average identity: {args["minid"]}%'))

        contig_ids = 'left_contig', 'right_contig'
        min_overlap_len = (args['min_overlap_len'] or
                           read_len//args['read_len_div'])
        visualize = args['visualize']
        seq = []
        start_inner_idx = 0
        if gaps[0][0] == 0:  # starting with a gap
            start_inner_idx = 1
            alignments = bam_util.fetch_alignments(*gaps[0])
            reads, read_positions = get_reads_with_positions(alignments)
            for read_id in reads:
                assembler = Assembler(
                        reads, read_positions, read_len,
                        {read_id: reads[read_id],
                         contig_ids[1]: contigs[0]},
                        (read_positions[read_id], gaps[0][1]+1),
                        min_overlap_len, visualize,
                        '%s_%d_%d' % (ref_id, *gaps[0]))
                part = assembler.assemble()
                if part:
                    seq.append(part)
                    break

        # Filling inner gaps
        for i in range(start_inner_idx, len(gaps)):
            # Checking if there is no right_contig contig (ending with a gap)
            if i+1-start_inner_idx > len(contigs)-1:
                break
            cutoff_len = 0 if not seq else len(contigs[i-start_inner_idx])
            alignments = bam_util.fetch_alignments(*gaps[i])
            reads, read_positions = get_reads_with_positions(alignments)
            assembler = Assembler(
                    reads, read_positions, read_len,
                    {contig_ids[0]: contigs[i-start_inner_idx],
                     contig_ids[1]: contigs[i+1-start_inner_idx]},
                    (gaps[i][0]-1, gaps[i][1]+1),
                    min_overlap_len, visualize,
                    '%s_%d_%d' % (ref_id, *gaps[i]))
            seq.append(assembler.assemble()[cutoff_len:])

        if gaps[-1][1] == ref_len-1:  # ending with a gap
            alignments = bam_util.fetch_alignments(*gaps[-1])
            reads, read_positions = get_reads_with_positions(alignments)
            for read_id in reversed(reads):
                assembler = Assembler(
                        reads, read_positions, read_len,
                        {contig_ids[0]: contigs[-1],
                         read_id: reads[read_id]},
                        (gaps[-1][0]-1, read_positions[read_id]),
                        min_overlap_len, visualize,
                        '%s_%d_%d' % (ref_id, *gaps[-1]))
                part = assembler.assemble()
                if part:
                    seq.append(part[:len(contigs[-1])])
                    break

        return ''.join(seq)


if __name__ == '__main__':
    parser = create_parser()
    args = vars(parser.parse_args())
    seq = main(args)
    print(seq)
    if args['output']:
        with open(args['output'], 'w', encoding='utf-8') as file:
            file.write('>output\n')
            file.write(seq)
