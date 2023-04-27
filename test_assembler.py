from argparse import ArgumentParser
import os

from Bio import SeqIO, Align

import assembler
import util


def create_parser():
    """Returns parser for command line arguments."""
    parser = ArgumentParser()
    parser.add_argument('max_tests_cnt', nargs='?', type=int,
                        help=('How many times at most to run assembler with '
                              'different data'))
    parser.add_argument('-r', '--read_len', type=int,
                        help='Only assemble sequences with reads of this '
                             'length')
    parser.add_argument('-mng', '--min_gap_len', type=int,
                        help='Only assemble sequences with gaps of length '
                             'greater or equal to this')
    parser.add_argument('-mxg', '--max_gap_len', type=int,
                        help='Only assemble sequences with gaps of length '
                             'shorter or equal to this')
    parser.add_argument('-d', '--read_len_divisor', type=int, default=3,
                        help=('How many times the minimum overlap length is '
                              'smaller than the read length '
                              '(used when min_overlap_len not given)'))
    parser.add_argument('-ol', '--min_overlap_len', type=int,
                        help='Minimum length of reads overlap')
    parser.add_argument('--test_data_dir', default='test_data/out_dir',
                        help='Directory containing simulated data')
    parser.add_argument('--summary', dest='summary', action='store_true',
                        help='Show only summary of testing')
    parser.set_defaults(summary=False)
    return parser


def test():
    parser = create_parser()
    args = parser.parse_args()

    description = ['Assembling simulated sequences with parameters:']
    if not args.read_len:
        description.append('read_len = any')
    else:
        description.append('read_len = %d' % args.read_len)
    if not args.max_gap_len:
        description.append('max_gap_len = any')
    else:
        description.append('max_gap_len = %d' % args.max_gap_len)
    if not args.min_overlap_len:
        description.append('min_overlap_len = read_len // %d' %
                           args.read_len_divisor)
    else:
        description.append('min_overlap_len = %d' % args.min_overlap_len)
    if args.max_tests_cnt:
        description.append('max_tests_cnt = %d' % args.max_tests_cnt)
    print('\n'.join(description), end='\n\n')

    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'

    tests_cnt = correct_cnt = wrong_cnt = failed_cnt = 0
    for data_dir in os.listdir(args.test_data_dir):
        if args.max_tests_cnt and tests_cnt >= args.max_tests_cnt:
            break

        # Expecting data_dir to be like 'sim_{read_len}_{gap_len}'
        read_len, gap_len = map(int, data_dir.split('_')[1:3])
        if (args.read_len and read_len != args.read_len
                or args.min_gap_len and gap_len < args.min_gap_len
                or args.max_gap_len and gap_len > args.max_gap_len):
            continue
        tests_cnt += 1

        if not args.summary:
            print(' | '.join((f'[Read length: {read_len}',
                              f'Gap length: {gap_len}',
                              f'Sequence length: {2*read_len + gap_len}]')))

        data_dir_path = os.path.join(os.path.dirname(__file__),
                                     args.test_data_dir, data_dir)
        target_file = f'{data_dir_path}/target.fa'
        target_seq = next(SeqIO.parse(target_file, 'fasta')).seq

        try:
            result_seq = assembler.main({
                'reads1': f'{data_dir_path}/dat1.fq',
                'reads2': f'{data_dir_path}/dat2.fq',
                'contig1': f'{data_dir_path}/cont1.fa',
                'contig2': f'{data_dir_path}/cont2.fa',
                'alignments_type': 'art',
                'alignments1': f'{data_dir_path}/dat1.aln',
                'alignments2': f'{data_dir_path}/dat2.aln',
                'read_len': read_len,
                'read_len_divisor': args.read_len_divisor,
                'min_overlap_len': args.min_overlap_len,
                'sam': args.sam
            })
            if result_seq == target_seq:
                correct_cnt += 1
                if not args.summary:
                    print('OK: resulting and target sequences match exactly')
            else:
                wrong_cnt += 1
                if not args.summary:
                    identity = util.compute_identity(
                        aligner, target_seq, result_seq)
                    print(
                        'Warning: resulting and target sequences differ',
                        f'- Identity: {identity}',
                        f'- Result: {result_seq}',
                        f'- Target: {target_seq}',
                        sep='\n'
                    )
        except assembler.WrongPathError:
            failed_cnt += 1
            if not args.summary:
                print('Error: Could not assemble the sequence')
        if not args.summary:
            print()

    summary = ['[SUMMARY]']
    summary.append(f'Correctly assembled: {correct_cnt}/{tests_cnt}')
    if wrong_cnt > 0:
        summary.append(f'Incorrectly assembled: {wrong_cnt}/{tests_cnt}')
    if failed_cnt > 0:
        summary.append(f'Could not assemble: {failed_cnt}/{tests_cnt}')
    print('\n'.join(summary))


if __name__ == '__main__':
    test()
