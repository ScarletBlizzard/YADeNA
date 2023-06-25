import itertools

from overlap_graph import OverlapGraph


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


def create_overlap_graph(reads, contigs, min_overlap_len,
                         visualize, graph_name):
    """
    Receives reads, contigs and minimum overlap length.
    Returns overlap graph created using hash table.
    """
    graph = OverlapGraph(name=graph_name, visualize=visualize)

    # Find overlaps between different reads, including between paired-end reads
    for r1, r2 in itertools.permutations(reads, 2):
        if r1 != r2:
            o_len = overlap_len(reads[r1], reads[r2], min_overlap_len)
            if o_len > 0:
                graph.add_child(r1, r2, o_len)

    c1, c2 = contigs.keys()
    for r in reads:
        # Find overlaps between suffix of left contig and prefixes of reads
        o_len = overlap_len(contigs[c1], reads[r], min_overlap_len)
        if o_len > 0:
            graph.add_child(c1, r, o_len)
        # Find overlaps between prefix of right contig and suffixes of reads
        o_len = overlap_len(reads[r], contigs[c2], min_overlap_len)
        if o_len > 0:
            graph.add_child(r, c2, o_len)

    return graph


class WrongPathError(Exception):
    """Raised when last read in path didn't match the read we want."""


def traverse(graph, reads, curr_read_id, last_read_id, visited=None):
    """
    Recursive function that traverses the overlap graph, thus assembling the
    target sequence. Uses greedy approach.

    Receives graph, current vertex (read), dict of reads, read that must be
    last in path and set of visited reads.
    Returns target sequence.
    """
    if visited is None:
        visited = set()
    visited.add(curr_read_id)
    while child := graph.get_next_child(curr_read_id):
        next_read_id, o_len = child

        if next_read_id in visited:
            continue
        try:
            res = traverse(graph, reads, next_read_id, last_read_id, visited)
            graph.mark(curr_read_id, next_read_id, o_len)
            return str(reads[curr_read_id]) + res[o_len:]
        except WrongPathError:
            continue
    if curr_read_id != last_read_id:
        raise WrongPathError('Could not assemble sequence')
    graph.mark(last_read_id)
    return str(reads[last_read_id])


def assemble(reads, contigs, min_overlap_len, visualize, graph_name):
    """Returns assembled pre-consensus sequence based on arguments."""
    if len(contigs) not in (1, 2):
        raise ValueError('Must have one or two contigs for assembly')

    contig_ids = list(contigs.keys())
    reads_with_contigs = {**reads, **contigs}
    min_min_overlap_len = min_overlap_len // 5
    while True and min_overlap_len > min_min_overlap_len:
        graph = create_overlap_graph(reads, contigs,
                                     min_overlap_len,
                                     visualize, graph_name)
        try:
            seq = traverse(graph, reads_with_contigs, *contig_ids)
            graph.visualize()
            return seq
        except WrongPathError:
            min_overlap_len -= 1
    return ''
