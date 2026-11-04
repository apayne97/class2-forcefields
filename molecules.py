"""Tiny test molecules (ethane, n-butane) with enumerated internal coordinates.

Coordinates are built from rough ideal geometries and then displaced randomly
(fixed seed) so that every bond, angle and torsion is away from its reference
value and all cross terms are non-zero. They are *not* optimised structures.
Lengths are in nm (OpenMM units).
"""

from dataclasses import dataclass, field
from itertools import combinations

import numpy as np

ANG = 0.1  # nm per angstrom


@dataclass
class Molecule:
    name: str
    types: list
    xyz: np.ndarray  # (n, 3), nm
    bonds: list  # list of (i, j), i < j
    neighbours: dict = field(init=False)

    def __post_init__(self):
        self.neighbours = {i: [] for i in range(len(self.types))}
        for i, j in self.bonds:
            self.neighbours[i].append(j)
            self.neighbours[j].append(i)

    def angles(self):
        """(a, vertex, c) with a < c."""
        out = []
        for b, nb in self.neighbours.items():
            for a, c in combinations(sorted(nb), 2):
                out.append((a, b, c))
        return out

    def torsions(self):
        """(a, b, c, d) for every b-c bond, each torsion once."""
        out = []
        for b, c in self.bonds:
            for a in self.neighbours[b]:
                for d in self.neighbours[c]:
                    if a != c and d != b and a != d:
                        out.append((a, b, c, d))
        return out

    def bond_pairs(self):
        """Pairs of bonds sharing an atom, as (a, b, c): bonds a-b and b-c."""
        return self.angles()

    def angle_pairs(self):
        """Pairs of angles with the same vertex and one shared neighbour.

        Returns (a, vertex, common, c): angles (a, vertex, common) and
        (c, vertex, common).
        """
        out = []
        for v, nb in self.neighbours.items():
            for common in nb:
                others = [x for x in nb if x != common]
                for a, c in combinations(sorted(others), 2):
                    out.append((a, v, common, c))
        return out


def _unit(v):
    return v / np.linalg.norm(v)


def _perp(v):
    t = np.array([1.0, 0, 0]) if abs(v[0]) < 0.9 else np.array([0, 1.0, 0])
    return _unit(np.cross(v, t))


def _methyl(c, to_neighbour, rch=1.111 * ANG, theta=111.0, phase=0.0):
    """Three H positions on carbon `c`, each at `theta` degrees from the C-C direction."""
    u = _unit(to_neighbour)
    p = _perp(u)
    q = np.cross(u, p)
    th = np.radians(theta)
    return [
        c + rch * (np.cos(th) * u + np.sin(th) * (np.cos(a) * p + np.sin(a) * q))
        for a in np.radians([phase, phase + 120, phase + 240])
    ]


def _methylene(c, d1, d2, rch=1.111 * ANG, alpha=54.75):
    """Two H positions on a carbon with two carbon neighbours in directions d1, d2."""
    u = _unit(-(_unit(d1) + _unit(d2)))
    n = _unit(np.cross(d1, d2))
    al = np.radians(alpha)
    return [c + rch * (np.cos(al) * u + s * np.sin(al) * n) for s in (1, -1)]


def _perturb(xyz, seed, scale=0.03 * ANG):
    rng = np.random.default_rng(seed)
    return xyz + rng.normal(scale=scale, size=xyz.shape)


def ethane(seed=0, perturb=True, staggered=True):
    cc = 1.535 * ANG
    c1, c2 = np.array([0, 0, -cc / 2]), np.array([0, 0, cc / 2])
    h1 = _methyl(c1, c2 - c1)
    h2 = _methyl(c2, c1 - c2, phase=0.0 if staggered else 60.0)
    xyz = np.array([c1, c2] + h1 + h2)
    types = ["C", "C", "H", "H", "H", "H", "H", "H"]
    bonds = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 5), (1, 6), (1, 7)]
    return Molecule("ethane", types, _perturb(xyz, seed) if perturb else xyz, bonds)


def butane(seed=0, perturb=True):
    """Trans n-butane, zig-zag backbone in the xy plane."""
    cc, th = 1.535 * ANG, np.radians(112.9)
    step = lambda s: cc * np.array([np.sin(th / 2), s * np.cos(th / 2), 0])  # noqa: E731
    c1 = np.zeros(3)
    c2 = c1 + step(+1)
    c3 = c2 + step(-1)
    c4 = c3 + step(+1)
    xyz = [c1, c2, c3, c4]
    types = ["C"] * 4
    bonds = [(0, 1), (1, 2), (2, 3)]
    for ci, nbrs in ((0, [1]), (3, [2])):
        for h in _methyl(xyz[ci], xyz[nbrs[0]] - xyz[ci]):
            bonds.append((ci, len(xyz)))
            xyz.append(h)
            types.append("H")
    for ci, (a, b) in ((1, (0, 2)), (2, (1, 3))):
        for h in _methylene(xyz[ci], xyz[a] - xyz[ci], xyz[b] - xyz[ci]):
            bonds.append((ci, len(xyz)))
            xyz.append(h)
            types.append("H")
    xyz = np.array(xyz)
    bonds = [(min(i, j), max(i, j)) for i, j in bonds]
    return Molecule("butane", types, _perturb(xyz, seed) if perturb else xyz, bonds)


MOLECULES = {"ethane": ethane, "butane": butane}
