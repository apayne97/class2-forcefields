"""Plain-numpy energies of every CFF93 valence term (the "hand-computed" reference).

Written independently of OpenMM, straight from Hwang et al. 1994 Eq. 1 and
Table 4. All energies in kcal/mol; lengths are converted from nm to angstrom.
Each sum counts every distinct bond/angle/torsion (or pair) once.
"""

import numpy as np

import cff93 as P

TERMS = [
    "bond", "angle", "torsion",
    "bond-bond", "bond-angle", "angle-angle",
    "bond-torsion", "angle-torsion", "angle-angle-torsion",
]


def length(x, i, j):
    return 10 * np.linalg.norm(x[i] - x[j])  # angstrom


def angle(x, a, b, c):
    u, v = x[a] - x[b], x[c] - x[b]
    return np.arccos(u @ v / np.linalg.norm(u) / np.linalg.norm(v))


def cos_dihedral(x, a, b, c, d):
    b1, b2, b3 = x[b] - x[a], x[c] - x[b], x[d] - x[c]
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    return n1 @ n2 / np.linalg.norm(n1) / np.linalg.norm(n2)


def dihedral(x, a, b, c, d):
    b1, b2, b3 = x[b] - x[a], x[c] - x[b], x[d] - x[c]
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    m = np.cross(n1, b2 / np.linalg.norm(b2))
    return np.arctan2(m @ n2, n1 @ n2)


def poly(d, k2, k3, k4):
    return k2 * d**2 + k3 * d**3 + k4 * d**4


def fourier(phi, k):
    return k[0] * np.cos(phi) + k[1] * np.cos(2 * phi) + k[2] * np.cos(3 * phi)


def types_of(mol, *idx):
    return P.canon(mol.types[i] for i in idx)


def bond_ref(mol, i, j):
    return P.BOND_T[(types_of(mol, i, j),)][0]


def angle_ref(mol, a, b, c):
    return np.radians(P.ANGLE_T[(types_of(mol, a, b, c),)][0])


def energies(mol):
    x = mol.xyz
    E = dict.fromkeys(TERMS, 0.0)
    for i, j in mol.bonds:
        b0, k2, k3, k4 = P.BOND_T[(types_of(mol, i, j),)]
        E["bond"] += poly(length(x, i, j) - b0, k2, k3, k4)
    for a, b, c in mol.angles():
        t0, k2, k3, k4 = P.ANGLE_T[(types_of(mol, a, b, c),)]
        E["angle"] += poly(angle(x, a, b, c) - np.radians(t0), k2, k3, k4)
    for a, b, c, d in mol.torsions():
        k = P.TORSION_T[(types_of(mol, a, b, c, d),)]
        phi = dihedral(x, a, b, c, d)
        E["torsion"] += sum(k[n] * (1 - np.cos((n + 1) * phi)) for n in range(3))

    # bond-bond: bonds a-b and b-c
    for a, b, c in mol.bond_pairs():
        key = tuple(sorted((types_of(mol, a, b), types_of(mol, b, c))))
        E["bond-bond"] += (
            P.BOND_BOND_T[key]
            * (length(x, a, b) - bond_ref(mol, a, b))
            * (length(x, b, c) - bond_ref(mol, b, c))
        )

    # bond-angle: each angle with each of its two bonds
    for a, b, c in mol.angles():
        akey = types_of(mol, a, b, c)
        dth = angle(x, a, b, c) - angle_ref(mol, a, b, c)
        for i, j in ((a, b), (b, c)):
            K = P.BOND_ANGLE_T[(types_of(mol, i, j), akey)]
            E["bond-angle"] += K * (length(x, i, j) - bond_ref(mol, i, j)) * dth

    # angle-angle: angles (a,v,common) and (c,v,common)
    for a, v, common, c in mol.angle_pairs():
        k1, k2 = types_of(mol, a, v, common), types_of(mol, c, v, common)
        K = P.ANGLE_ANGLE_T[(tuple(sorted((k1, k2))), types_of(mol, v, common))]
        E["angle-angle"] += (
            K
            * (angle(x, a, v, common) - angle_ref(mol, a, v, common))
            * (angle(x, c, v, common) - angle_ref(mol, c, v, common))
        )

    for a, b, c, d in mol.torsions():
        tkey = types_of(mol, a, b, c, d)
        phi = dihedral(x, a, b, c, d)
        # bond-torsion: central bond and both terminal bonds
        K = P.BOND_TORSION_CENTRAL_T[(types_of(mol, b, c), tkey)]
        E["bond-torsion"] += (length(x, b, c) - bond_ref(mol, b, c)) * fourier(phi, K)
        for i, j in ((a, b), (c, d)):
            K = P.BOND_TORSION_TERMINAL_T[(types_of(mol, i, j), tkey)]
            E["bond-torsion"] += (length(x, i, j) - bond_ref(mol, i, j)) * fourier(phi, K)
        # angle-torsion: both angles contained in the torsion
        for p, q, r in ((a, b, c), (b, c, d)):
            K = P.ANGLE_TORSION_T[(types_of(mol, p, q, r), tkey)]
            E["angle-torsion"] += (angle(x, p, q, r) - angle_ref(mol, p, q, r)) * fourier(phi, K)
        # angle-angle-torsion
        k1, k2 = types_of(mol, a, b, c), types_of(mol, b, c, d)
        K = P.ANGLE_ANGLE_TORSION_T[(tuple(sorted((k1, k2))), tkey)]
        E["angle-angle-torsion"] += (
            K
            * (angle(x, a, b, c) - angle_ref(mol, a, b, c))
            * (angle(x, b, c, d) - angle_ref(mol, b, c, d))
            * np.cos(phi)
        )
    return E
