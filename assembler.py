import itertools

from overlap_graph import OverlapGraph, GraphVisualizer


class WrongPathError(Exception):
    """Raised when last read in path didn't match the read we want."""


class Assembler:
    def __init__(self, reads: dict[str, str], read_len: int,
                 contigs: dict[str, str], min_overlap_len: int,
                 visualize: bool, graph_name: str, positions=None):
        if len(contigs) not in (1, 2):
            raise ValueError('Must have one or two contigs for assembly')

        self._reads = reads
        self._read_len = read_len
        self._min_overlap_len = min_overlap_len

        self._contig_ids = tuple(contigs.keys())
        self._contigs = tuple(contigs.values())

        self._visualizer = GraphVisualizer(
                graph_name, positions, 'pos: ') if visualize else None

        self._graph_name = graph_name  # For debug

    def assemble(self):
        while True and self._min_overlap_len > 0:
            self._create_overlap_graph()
            if self._visualizer:
                for cid in self._contig_ids:
                    self._visualizer.mark_node(cid)
            try:
                seq = self._traverse_overlap_graph(self._contig_ids[0])
                if self._visualizer:
                    self._visualizer.render()
                return seq
            except WrongPathError:
                self._min_overlap_len -= 1
        if self._visualizer:
            self._visualizer.render()
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
        if self._visualizer:
            self._visualizer.empty()

        # Find overlaps of reads
        for r1, r2 in itertools.permutations(self._reads, 2):
            o_len = self._overlap_len(self._reads[r1], self._reads[r2])
            if o_len > 0:
                self._graph.add_child(r1, r2, o_len)
                if self._visualizer:
                    self._visualizer.add_child(r1, r2, str(o_len))

        # Find overlaps between suffix of left contig and prefixes of reads
        for r in self._reads:
            o_len = self._overlap_len(self._contigs[0], self._reads[r])
            if o_len > 0:
                self._graph.add_child(self._contig_ids[0], r, o_len)
                if self._visualizer:
                    self._visualizer.add_child(
                            self._contig_ids[0], r, str(o_len))

        # Find overlaps between prefix of right contig and suffixes of reads
        for r in self._reads:
            o_len = self._overlap_len(self._reads[r], self._contigs[1])
            if o_len > 0:
                self._graph.add_child(r, self._contig_ids[1], o_len)
                if self._visualizer:
                    self._visualizer.add_child(
                            r, self._contig_ids[1], str(o_len))

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
                if self._visualizer:
                    self._visualizer.mark_connection(
                            curr_read_id, next_read_id, str(o_len))
                if curr_read_id != self._contig_ids[0]:
                    return self._reads[curr_read_id] + res
                else:
                    return res
            except WrongPathError:
                continue
        if curr_read_id != self._contig_ids[1]:
            raise WrongPathError('Could not assemble sequence')
        return ''
