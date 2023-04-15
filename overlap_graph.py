from collections import defaultdict
import heapq


class OverlapGraph:
    def __init__(self):
        self.__graph = defaultdict(list)

    def add_child(self, parent_id: str, child_id: str, overlap_len: int):
        heapq.heappush(self.__graph[parent_id], (-overlap_len, child_id))

    def get_next_child(self, parent_id):
        try:
            child = heapq.heappop(self.__graph[parent_id])
            return (child[1], -child[0])
        except IndexError:
            return None
