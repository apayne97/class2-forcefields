"""CFF93 hydrocarbon parameters, transcribed from Hwang, Stockfisch & Hagler,
J. Am. Chem. Soc. 116, 2515 (1994), Table 4 (pp. 2520-2521).

Units are those printed in the paper: kcal/mol, angstrom, radian, degree.
These are the *published CFF93 alkane values*; they are used here to build
small worked examples, not as a production parameter set. Nonbonded terms
(Table 4 D and E) are not included.

Keys are strings written as in the paper ("H-C-C/H-C-C-H"). `canon` makes them
direction-independent so that "H-C-C" and "C-C-H" are the same angle.
"""

# Table 4 A: E = 2K (b-b0)^2 + 3K (b-b0)^3 + 4K (b-b0)^4
# value: (b0 [A], 2K [kcal/mol/A^2], 3K [kcal/mol/A^3], 4K [kcal/mol/A^4])
BOND = {
    "H-C": (1.111, 346.4, -706.6, 863.2),
    "C-C": (1.535, 299.4, -515.8, 667.1),
}

# Table 4 B: same form in theta. value: (theta0 [deg], 2K, 3K, 4K [kcal/mol/rad^n])
ANGLE = {
    "H-C-H": (107.7, 41.7, -7.2, -8.8),
    "H-C-C": (110.8, 42.7, -8.8, -9.1),
    "C-C-C": (112.9, 42.3, -9.8, -9.2),
}

# Table 4 C: E = 1K(1-cos phi) + 2K(1-cos 2phi) + 3K(1-cos 3phi) [kcal/mol]
TORSION = {
    "H-C-C-H": (-0.968, 0.010, -0.150),
    "H-C-C-C": (-0.968, 0.010, -0.150),
    "C-C-C-C": (-0.968, 0.010, -0.150),
}

# Table 4 F: E = Kbb' (b-b0)(b'-b0')  [kcal/mol/A^2]
BOND_BOND = {
    "H-C/H-C": 10.24,
    "H-C/C-C": 9.84,
    "C-C/C-C": 8.75,
}

# Table 4 G: E = Kbtheta (b-b0)(theta-theta0)  [kcal/mol/A/rad]
BOND_ANGLE = {
    "H-C/H-C-H": 19.7,
    "H-C/H-C-C": 12.1,
    "C-C/H-C-C": 30.0,
    "C-C/C-C-C": 16.0,
}

# Table 4 H: E = Ktt' (theta-theta0)(theta'-theta0')  [kcal/mol/rad^2]
# key: "angle/angle|common bond" (the two angles share a vertex and this bond)
ANGLE_ANGLE = {
    "H-C-H/H-C-H|H-C": 0.81,
    "H-C-H/H-C-C|H-C": 1.36,
    "H-C-C/H-C-C|H-C": 2.71,
    "H-C-C/H-C-C|C-C": -1.90,
    "H-C-C/C-C-C|C-C": -4.02,
    "C-C-C/C-C-C|C-C": -7.81,
}

# Table 4 I: E = Kphitt' (theta-theta0)(theta'-theta0') cos(phi)  [kcal/mol/rad^2]
ANGLE_ANGLE_TORSION = {
    "H-C-C/H-C-C/H-C-C-H": -12.9,
    "H-C-C/C-C-C/H-C-C-C": -16.7,
    "C-C-C/C-C-C/C-C-C-C": -31.0,
}

# Table 4 J: E = (b-b0)[1K cos phi + 2K cos 2phi + 3K cos 3phi], b = central bond
# [kcal/mol/A]
BOND_TORSION_CENTRAL = {
    "C-C/H-C-C-H": (-45.4, -0.68, -0.78),
    "C-C/H-C-C-C": (-42.7, -2.76, -0.34),
    "C-C/C-C-C-C": (-41.1, -4.91, 0.27),
}

# Table 4 K: same form, b' = a terminal bond of the torsion [kcal/mol/A]
BOND_TORSION_TERMINAL = {
    "H-C/H-C-C-H": (0.68, 0.41, 0.07),
    "H-C/H-C-C-C": (0.68, 0.25, 0.35),
    "C-C/H-C-C-C": (1.45, 0.51, -0.30),
    "C-C/C-C-C-C": (2.00, 0.67, 0.31),
}

# Table 4 L: E = (theta-theta0)[1K cos phi + 2K cos 2phi + 3K cos 3phi]
# [kcal/mol/rad]; theta is either angle that is part of the torsion
ANGLE_TORSION = {
    "H-C-C/H-C-C-H": (-1.41, 0.64, -0.39),
    "H-C-C/H-C-C-C": (0.23, 0.61, -0.14),
    "C-C-C/H-C-C-C": (-1.64, 0.09, -0.26),
    "C-C-C/C-C-C-C": (-0.03, -0.01, 0.00),
}


def canon(chain):
    """Direction-independent key for a chain of atom types, e.g. ('C','C','H')."""
    chain = tuple(chain)
    return min(chain, chain[::-1])


def parse(key):
    """'H-C-C/H-C-C-H' -> (('C','C','H'), ('C','C','H','H')...) canonical chains."""
    return tuple(canon(part.split("-")) for part in key.split("/"))


def _table(raw):
    return {parse(k): v for k, v in raw.items()}


BOND_T = _table(BOND)
ANGLE_T = _table(ANGLE)
TORSION_T = _table(TORSION)
BOND_BOND_T = {tuple(sorted(parse(k))): v for k, v in BOND_BOND.items()}
BOND_ANGLE_T = _table(BOND_ANGLE)
ANGLE_TORSION_T = _table(ANGLE_TORSION)
BOND_TORSION_CENTRAL_T = _table(BOND_TORSION_CENTRAL)
BOND_TORSION_TERMINAL_T = _table(BOND_TORSION_TERMINAL)
ANGLE_ANGLE_TORSION_T = {
    (tuple(sorted(parse(k.rsplit("/", 1)[0]))), parse(k.rsplit("/", 1)[1])[0]): v
    for k, v in ANGLE_ANGLE_TORSION.items()
}


def _aa_table(raw):
    out = {}
    for k, v in raw.items():
        angles, bond = k.split("|")
        a1, a2 = angles.split("/")
        key = (tuple(sorted((canon(a1.split("-")), canon(a2.split("-"))))), canon(bond.split("-")))
        out[key] = v
    return out


ANGLE_ANGLE_T = _aa_table(ANGLE_ANGLE)
