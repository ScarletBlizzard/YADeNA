from argparse import ArgumentParser

from Bio import SeqIO

import assembler
from util import BAMUtil


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


def find_gaps(bam_util, ref_len, read_len, min_avg_identity):
    gaps = []
    window_len = 2 * read_len
    start, step = 0, read_len
    while start + step < ref_len:
        end = min(start+window_len-1, ref_len-1)
        avg_identity = bam_util.compute_avg_identity(start, end)
        if avg_identity < min_avg_identity:
            if not gaps or start-gaps[-1][1] > step+1:
                gaps.append((start, end))
            else:
                prev_start, _ = gaps.pop()
                gaps.append((prev_start, end))
        start += step
    return gaps


def get_contigs(bam_util, gaps, ref_len):
    if gaps:
        contigs = []
        if gaps[0][0] != 0:  # starting with a contig
            contigs.append(bam_util.consensus(0, gaps[0][0]-1))
        for i in range(1, len(gaps)):
            contigs.append(bam_util.consensus(gaps[i-1][1]+1, gaps[i][0]-1))
        if gaps[-1][1] != ref_len-1:  # ending with a contig
            contigs.append(bam_util.consensus(gaps[-1][1]+1, ref_len-1))
        return contigs
    return [bam_util.consensus(0, ref_len-1)]


class IdentityError(ValueError):
    pass


def main(args):
    ref_record = next(SeqIO.parse(args['reference'], 'fasta'))
    ref_id, ref_len = ref_record.name, len(ref_record.seq)
    read_len = args['read_len']
    with BAMUtil(args['reference'], ref_id,
                 [args['reads1'], args['reads2']]) as bam_util:
        gaps = find_gaps(
                bam_util, ref_len, read_len, args['minid']/100)
        contigs = get_contigs(bam_util, gaps, ref_len)
        if not gaps:
            return contigs[0]
        if not contigs:
            raise IdentityError(
                    ('Cannot assemble sequence with given '
                     f'minimum average identity: {args["minid"]}%'))

        lcontig_id, rcontig_id = 'left_contig', 'right_contig'
        min_overlap_len = (args['min_overlap_len'] or
                           read_len//args['read_len_div'])
        seq = []
        start_inner_idx = 0
        if gaps[0][0] == 0:  # starting with a gap
            start_inner_idx = 1
            reads, read_positions = bam_util.fetch_reads(*gaps[0])
            read_positions[rcontig_id] = gaps[0][1] + 1
            for read_id in reads:
                part = assembler.assemble(
                        reads, read_positions, read_len,
                        {read_id: reads[read_id],
                         rcontig_id: contigs[0]},
                        min_overlap_len, args['visualize'],
                        graph_name='Gap_%d_%d' % gaps[0])
                if part:
                    seq.append(part)
                    break

        # Filling inner gaps
        for i in range(start_inner_idx, len(gaps)):
            # Checking if there is no right_contig contig (ending with a gap)
            if i+1-start_inner_idx > len(contigs)-1:
                break
            cutoff_len = 0 if not seq else len(contigs[i-start_inner_idx])
            reads, read_positions = bam_util.fetch_reads(*gaps[i])
            read_positions[lcontig_id] = gaps[i][0] - 1
            read_positions[rcontig_id] = gaps[i][1] + 1
            seq.append(assembler.assemble(
                    reads, read_positions, read_len,
                    {lcontig_id: contigs[i-start_inner_idx],
                     rcontig_id: contigs[i+1-start_inner_idx]},
                    min_overlap_len, args['visualize'],
                    graph_name='Gap_%d_%d' % gaps[i]
                    )[cutoff_len:])

        if gaps[-1][1] == ref_len-1:  # ending with a gap
            reads, read_positions = bam_util.fetch_reads(*gaps[-1])
            read_positions[lcontig_id] = gaps[-1][0] - 1
            for read_id in reversed(reads):
                part = assembler.assemble(
                    reads, read_positions, read_len,
                    {lcontig_id: contigs[-1],
                     read_id: reads[read_id]},
                    min_overlap_len, args['visualize'],
                    graph_name='Gap_%d_%d' % gaps[-1])
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
            file.write('>output\n')  # TODO: Change this string
            file.write(seq)
