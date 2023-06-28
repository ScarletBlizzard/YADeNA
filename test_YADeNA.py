from argparse import ArgumentParser
import os
import tempfile
import traceback

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import util
import YADeNA


def create_parser():
    """Returns parser for command line arguments."""
    parser = ArgumentParser()
    parser.add_argument('-i', '--minid', required=True, type=float,
                        help=('Minimum average percent-identity of a '
                              'genome region to consider it a gap'))
    parser.add_argument('max_tests_cnt', nargs='?', type=int,
                        help=('How many times at most to run assembler with '
                              'different data'))
    parser.add_argument('-l', '--read_len', type=int,
                        help='Only assemble sequences with reads of this '
                             'length')
    parser.add_argument('-mn', '--min_seq_len', type=int,
                        help='Only assemble sequences with gaps of length '
                             'greater or equal to this')
    parser.add_argument('-mx', '--max_seq_len', type=int,
                        help='Only assemble sequences with gaps of length '
                             'shorter or equal to this')
    parser.add_argument('-d', '--read_len_div', type=int, default=3,
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
    """Function for testing assembler on simulated data."""
    parser = create_parser()
    args = parser.parse_args()

    description = ['Assembling simulated sequences with parameters:']
    if not args.read_len:
        description.append('read_len = any')
    else:
        description.append('read_len = %d' % args.read_len)
    if not args.max_seq_len:
        description.append('max_seq_len = any')
    else:
        description.append('max_seq_len = %d' % args.max_seq_len)
    if not args.min_overlap_len:
        description.append('min_overlap_len = read_len // %d' %
                           args.read_len_div)
    else:
        description.append('min_overlap_len = %d' % args.min_overlap_len)
    if args.max_tests_cnt:
        description.append('max_tests_cnt = %d' % args.max_tests_cnt)
    print('\n'.join(description), end='\n\n')

    tests_cnt = assembled_cnt = failed_cnt = identity_sum = 0
    for data_dir in os.listdir(args.test_data_dir):
        if args.max_tests_cnt and tests_cnt >= args.max_tests_cnt:
            break

        # Expecting data_dir to be like 'sim_{read_len}_{seq_len}'
        read_len, seq_len = map(int, data_dir.split('_')[1:3])
        if (args.read_len and read_len != args.read_len
                or args.min_seq_len and seq_len < args.min_seq_len
                or args.max_seq_len and seq_len > args.max_seq_len):
            continue
        tests_cnt += 1

        if not args.summary:
            print(' | '.join((f'[Read length: {read_len}',
                              f'Sequence length: {seq_len}]')))

        data_dir_path = os.path.join(os.path.dirname(__file__),
                                     args.test_data_dir, data_dir)
        target_fname = f'{data_dir_path}/target.fa'

        try:
            result_seq = YADeNA.main({
                'reference': target_fname,
                'minid': args.minid,
                'reads1': f'{data_dir_path}/dat1.fq',
                'reads2': f'{data_dir_path}/dat2.fq',
                'read_len': read_len,
                'read_len_div': args.read_len_div,
                'min_overlap_len': args.min_overlap_len,
                'visualize': False
            })
            assembled_cnt += 1
            result_fname = None
            with tempfile.NamedTemporaryFile(mode='w', delete=False) as r_file:
                SeqIO.write(SeqRecord(Seq(result_seq), id='result'),
                            r_file, 'fasta')
                result_fname = r_file.name
            identity = util.compute_identity_of_unaligned(
                    target_fname, result_fname)
            os.remove(result_fname)
            print(f'Identity: {identity}')
            identity_sum += identity
        except Exception:
            failed_cnt += 1
            if not args.summary:
                print(traceback.format_exc())
        if not args.summary:
            print()

    summary = ['[SUMMARY]']
    summary.append(f'Assembled: {assembled_cnt}/{tests_cnt}')
    if failed_cnt > 0:
        summary.append(f'Could not assemble: {failed_cnt}/{tests_cnt}')
    summary.append(f'Average identity: {identity_sum/assembled_cnt}')
    print('\n'.join(summary))


if __name__ == '__main__':
    test()
