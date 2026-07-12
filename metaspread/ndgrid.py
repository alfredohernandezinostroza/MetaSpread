"""A minimal N-dimensional multi-occupancy grid.

mesa 2.1.2 ships only 2D spaces, so the optional 3D simulation path uses this
lightweight grid instead. It implements just the subset of the mesa
``MultiGrid`` API that MetaSpread actually calls — ``place_agent``,
``move_agent``, ``remove_agent``, ``get_neighborhood`` (von Neumann and Moore),
``get_cell_list_contents`` and ``out_of_bounds`` — for an arbitrary number of
non-toroidal dimensions. The 2D path keeps using mesa unchanged, so this grid is
only exercised when ``space_dimensions == 3``.
"""
from collections import defaultdict
from itertools import product


class NDGrid:
    def __init__(self, dimensions):
        self.dimensions = tuple(int(d) for d in dimensions)
        self.ndim = len(self.dimensions)
        self.width = self.dimensions[0]
        self.height = self.dimensions[1]
        self.depth = self.dimensions[2] if self.ndim > 2 else None
        self._grid = defaultdict(list)

    def out_of_bounds(self, pos):
        return any(c < 0 or c >= d for c, d in zip(pos, self.dimensions))

    def place_agent(self, agent, pos):
        pos = tuple(pos)
        self._grid[pos].append(agent)
        agent.pos = pos

    def remove_agent(self, agent):
        pos = tuple(agent.pos)
        self._grid[pos].remove(agent)
        if not self._grid[pos]:
            del self._grid[pos]
        agent.pos = None

    def move_agent(self, agent, pos):
        self.remove_agent(agent)
        self.place_agent(agent, pos)

    def get_cell_list_contents(self, cell_list):
        contents = []
        for pos in cell_list:
            contents.extend(self._grid.get(tuple(pos), ()))
        return contents

    def get_neighborhood(self, pos, moore=False, include_center=False):
        pos = tuple(pos)
        neighborhood = []
        if moore:
            for delta in product((-1, 0, 1), repeat=self.ndim):
                if not include_center and all(d == 0 for d in delta):
                    continue
                candidate = tuple(c + d for c, d in zip(pos, delta))
                if not self.out_of_bounds(candidate):
                    neighborhood.append(candidate)
        else:
            if include_center:
                neighborhood.append(pos)
            for axis in range(self.ndim):
                for step in (-1, 1):
                    candidate = list(pos)
                    candidate[axis] += step
                    candidate = tuple(candidate)
                    if not self.out_of_bounds(candidate):
                        neighborhood.append(candidate)
        return neighborhood
