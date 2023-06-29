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
    match_count = 0
    i, j, d = 0, 0, 0
    for op, ln in cigartuples:
        if op == 0:  # M
            for _ in range(ln):
                if read[i-d] == seq_part[j+d]:
                    match_count += 1
                i += 1
                j += 1
        else:
            i += ln
            if op == 2:  # D
                d += ln
    return match_count / i


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


def get_diff_segments(segments, end):
    diff = []
    if segments:
        if segments[0][0] != 0:
            diff.append((0, segments[0][0]-1))
        for i in range(1, len(segments)):
            diff.append((segments[i-1][1]+1, segments[i][0]-1))
        if segments[-1][1] != end:
            diff.append((segments[-1][1]+1, end))
    return diff


def get_contigs(contig_segments, consensus_seq):
    contigs = []
    for seg in contig_segments:
        contigs.append(consensus_seq[seg[0]:seg[1]+1])
    return contigs


def get_reads(alignments: Alignments):
    reads = {read_id: aln.query_sequence
             for read_id, aln in alignments.items()}
    return reads


def get_read_positions(alignments: Alignments):
    positions = {read_id: aln.reference_start
                 for read_id, aln in alignments.items()}
    return positions


class CommandArgumentsError(Exception):
    """Exception raised for errors in the command line args."""


def main(args):
    ref_record = next(SeqIO.parse(args['reference'], 'fasta'))
    ref_id, ref_len = ref_record.name, len(ref_record.seq)
    read_len = args['read_len']
    seq = []

    with BAMUtil(args['reference'], [args['reads1'],
                 args['reads2']]) as bam_util:
        consensus_seq = bam_util.consensus()
        gaps = find_gaps(
                bam_util, consensus_seq, read_len, args['minid']/100)
        if not gaps:
            return consensus_seq

        contig_segments = get_diff_segments(gaps, len(consensus_seq)-1)
        contigs = get_contigs(contig_segments, consensus_seq)
        if not contigs:
            raise CommandArgumentsError(
                    ('Cannot assemble sequence with given '
                     f'minimum average identity: {args["minid"]}%'))

        contig_ids = 'left_contig', 'right_contig'
        min_overlap_len = (args['min_overlap_len'] or
                           read_len//args['read_len_div'])
        visualize = args['visualize']
        start_inner_idx = 0
        if gaps[0][0] == 0:  # starting with a gap
            start_inner_idx = 1
            alignments = bam_util.fetch_alignments(*gaps[0])
            reads = get_reads(alignments)
            positions = get_read_positions(alignments)
            positions[contig_ids[1]] = contig_segments[0][0]
            for read_id in reads:
                assembler = Assembler(
                        reads, read_len,
                        {read_id: reads[read_id],
                         contig_ids[1]: contigs[0]},
                        min_overlap_len, visualize,
                        '%s_%d_%d' % (ref_id, *gaps[0]),
                        positions)
                part = assembler.assemble()
                if part:
                    seq.append(part)
                    break

        # Filling inner gaps
        for i in range(start_inner_idx, len(gaps)):
            # Checking if there is no right_contig contig (ending with a gap)
            contig_idx = i-start_inner_idx
            if contig_idx+1 > len(contigs)-1:
                break
            alignments = bam_util.fetch_alignments(*gaps[i])
            positions = get_read_positions(alignments)
            positions[contig_ids[0]] = contig_segments[contig_idx][0]
            positions[contig_ids[1]] = contig_segments[contig_idx+1][0]
            assembler = Assembler(
                    get_reads(alignments), read_len,
                    {contig_ids[0]: contigs[contig_idx],
                     contig_ids[1]: contigs[contig_idx+1]},
                    min_overlap_len, visualize,
                    '%s_%d_%d' % (ref_id, *gaps[i]),
                    positions)
            seq.append(contigs[contig_idx])
            seq.append(assembler.assemble())

        seq.append(contigs[-1])

        if gaps[-1][1] == ref_len-1:  # ending with a gap
            alignments = bam_util.fetch_alignments(*gaps[-1])
            reads = get_reads(alignments)
            positions = get_read_positions(alignments)
            positions[contig_ids[0]] = contig_segments[-1][0]
            for read_id in reversed(reads):
                assembler = Assembler(
                        reads, read_len,
                        {contig_ids[0]: contigs[-1],
                         read_id: reads[read_id]},
                        min_overlap_len, visualize,
                        '%s_%d_%d' % (ref_id, *gaps[-1]),
                        positions)
                part = assembler.assemble()
                if part:
                    seq.append(part)
                    break

    with BAMUtil.from_ref_str(
            ''.join(seq), 'output',
            [args['reads1'], args['reads2']]) as bam_util:
        bam_util.save_sam()
        return bam_util.consensus()


if __name__ == '__main__':
    parser = create_parser()
    args = vars(parser.parse_args())
    try:
        seq = main(args)
        print(seq)
        if args['output']:
            with open(args['output'], 'w', encoding='utf-8') as file:
                file.write('>output\n')
                file.write(seq)
    except CommandArgumentsError as e:
        print(e)
