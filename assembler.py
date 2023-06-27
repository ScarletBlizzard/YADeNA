import itertools

from frozendict import frozendict

from overlap_graph import OverlapGraph, GraphVisualizer


class WrongPathError(Exception):
    """Raised when last read in path didn't match the read we want."""


class Assembler:
    def __init__(self, sorted_reads: frozendict[str, str],
                 read_positions: dict[str, int], read_len: int,
                 contigs: dict[str, str], contig_positions: tuple[int, ...],
                 min_overlap_len: int, visualize: bool, graph_name: str):
        if len(contigs) not in (1, 2):
            raise ValueError('Must have one or two contigs for assembly')

        self._sorted_reads = sorted_reads
        self._read_positions = read_positions
        self._read_len = read_len
        self._min_overlap_len = min_overlap_len

        self._contig_ids = tuple(contigs.keys())
        self._contigs = tuple(contigs.values())
        self._contig_positions = contig_positions

        for contig_id, pos in zip(self._contig_ids, contig_positions):
            self._read_positions[contig_id] = pos

        self._graph_visualizer = GraphVisualizer(
                    graph_name, read_positions, 'pos: '
                    ) if visualize else None

    def assemble(self):
        min_min_overlap_len = self._min_overlap_len // 5
        while True and self._min_overlap_len > min_min_overlap_len:
            self._create_overlap_graph()
            try:
                seq = self._traverse_overlap_graph(self._contig_ids[0])
                if self._graph_visualizer:
                    self._graph_visualizer.render()
                return seq
            except WrongPathError:
                self._min_overlap_len -= 1
        return ''

    def _overlap_len(self, s1: str, s2: str):
        """
        Receives strings s1 and s2.
        Returns the length of the longest suffix of s1 that matches
        the prefix of s2.
        """
        start = 0
        while True:
            start = s1.find(s2[:self._min_overlap_len], start)
            if start == -1:
                return 0
            if s2.startswith(s1[start:]):
                return len(s1)-start
            start += 1

    def _create_overlap_graph(self):
        self._graph = OverlapGraph()
        max_dist = self._read_len - self._min_overlap_len

        # Find overlaps of reads
        for r1, r2 in itertools.combinations(self._sorted_reads, 2):
            if self._read_positions[r2]-self._read_positions[r1] > max_dist:
                continue
            o_len = self._overlap_len(self._sorted_reads[r1],
                                      self._sorted_reads[r2])
            if o_len > 0:
                self._graph.add_child(r1, r2, o_len)
                if self._graph_visualizer:
                    self._graph_visualizer.add_child(r1, r2, str(o_len))

        # Find overlaps between suffix of left contig and prefixes of reads
        for r in self._sorted_reads:
            if self._read_positions[r]-self._contig_positions[0] > max_dist:
                break
            o_len = self._overlap_len(self._contigs[0], self._sorted_reads[r])
            if o_len > 0:
                self._graph.add_child(self._contig_ids[0], r, o_len)
                if self._graph_visualizer:
                    self._graph_visualizer.add_child(self._contig_ids[0],
                                                     r, str(o_len))

        # Find overlaps between prefix of right contig and suffixes of reads
        for r in reversed(self._sorted_reads):
            if self._contig_positions[1]-self._read_positions[r] > max_dist:
                break
            o_len = self._overlap_len(self._sorted_reads[r], self._contigs[1])
            if o_len > 0:
                self._graph.add_child(r, self._contig_ids[1], o_len)
                if self._graph_visualizer:
                    self._graph_visualizer.add_child(r, self._contig_ids[1],
                                                     str(o_len))

    def _traverse_overlap_graph(self, curr_read_id: str, visited=None):
        if visited is None:
            visited = set()
        visited.add(curr_read_id)
        while child := self._graph.get_next_child(curr_read_id):
            next_read_id, o_len = child

            if next_read_id in visited:
                continue
            try:
                res = self._traverse_overlap_graph(
                    next_read_id, visited)[o_len:]
                if self._graph_visualizer:
                    self._graph_visualizer.mark(
                            curr_read_id, next_read_id, str(o_len))
                try:
                    return self._sorted_reads[curr_read_id] + res
                except KeyError:  # occurs only if id is that of left contig
                    return self._contigs[0] + res
            except WrongPathError:
                continue
        if curr_read_id != self._contig_ids[1]:
            raise WrongPathError('Could not assemble sequence')
        if self._graph_visualizer:
            self._graph_visualizer.mark(self._contig_ids[1])
        return str(self._contigs[1])
