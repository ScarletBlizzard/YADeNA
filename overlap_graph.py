from collections import defaultdict
import heapq

import graphviz


class OverlapGraph:
    def __init__(self, graph=None):
        self._graph = graph or defaultdict(list)

    def add_child(self, parent_id: str, child_id: str, overlap_len: int):
        heapq.heappush(self._graph[parent_id], (-overlap_len, child_id))

    def get_next_child(self, parent_id: str):
        try:
            child = heapq.heappop(self._graph[parent_id])
            return (child[1], -child[0])
        except IndexError:
            return None

    def __repr__(self):
        return f'OverlapGraph({self._graph})'


class GraphVisualizer:
    def __init__(self, name, labels=None, label_prefix=''):
        self._name = name + '_graph.gv'
        self._labels = labels
        self._label_prefix = label_prefix
        self._font = 'helvetica'
        self._mark_color = 'red'
        self._graph_renderer = graphviz.Digraph(self._name, strict=True)

    def _get_label(self, node_id: str, color='black'):
        if not self._labels:
            return ''
        return (f'<<font color="{color}">'
                f'{self._label_prefix}{self._labels[node_id]}'
                '</font>>')

    def _add_node(self, node_id: str):
        self._graph_renderer.node(node_id, node_id, fontname=self._font,
                                  xlabel=self._get_label(node_id))

    def _add_edge(self, from_id: str, to_id: str, label: str):
        self._graph_renderer.edge(from_id, to_id, label,
                                  fontname=self._font)

    def add_child(self, parent_id: str, child_id: str, label: str):
        self._add_node(parent_id)
        self._add_node(child_id)
        self._add_edge(parent_id, child_id, label)

    def mark_node(self, node_id: str):
        self._graph_renderer.node(
                node_id, node_id, style='filled',
                color=self._mark_color, fontcolor='white',
                fontname=f'{self._font}-bold',
                xlabel=self._get_label(node_id, self._mark_color))

    def mark_connection(self, parent_id: str, child_id=None, label=None):
        self.mark_node(parent_id)
        if child_id and label:
            self._graph_renderer.edge(
                    parent_id, child_id, label,
                    color=self._mark_color,
                    fontcolor=self._mark_color)

    def render(self):
        self._graph_renderer.render(self._name, view=False)

    def empty(self):
        self._graph_renderer = graphviz.Digraph(self._name, strict=True)
