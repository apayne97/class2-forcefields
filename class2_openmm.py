"""OpenMM implementation of the CFF93 valence terms (Hwang et al. 1994, Eq. 1).

One OpenMM force per term family, each in its own force group so the energy of
each term can be read separately:

    term                  OpenMM class                 expression variables
    --------------------  ---------------------------  ---------------------
    bond (2K,3K,4K)       CustomBondForce              r
    angle (2K,3K,4K)      CustomAngleForce             theta
    torsion (1K,2K,3K)    CustomTorsionForce           theta  (the dihedral)
    bond-bond             CustomCompoundBondForce (3)  distance(p1,p2), distance(p2,p3)
    bond-angle            CustomCompoundBondForce (3)  distance, angle
    angle-angle           CustomCompoundBondForce (4)  angle(p1,p2,p3), angle(p4,p2,p3)
    bond-torsion          CustomCompoundBondForce (4)  distance, dihedral
    angle-torsion         CustomCompoundBondForce (4)  angle, dihedral
    angle-angle-torsion   CustomCompoundBondForce (4)  angle, angle, dihedral

Parameters are the published CFF93 alkane values (cff93.py), converted from
kcal/mol, angstrom, rad to OpenMM units (kJ/mol, nm, rad).
"""

import click
import numpy as np
import openmm as mm

import class2_reference as ref
import cff93 as P
from molecules import MOLECULES

KJ = 4.184  # kJ per kcal
NM = 0.1  # nm per angstrom
TERM_GROUP = {name: i for i, name in enumerate(ref.TERMS)}


def _ty(mol, *idx):
    return P.canon(mol.types[i] for i in idx)


def _b0(mol, i, j):  # nm
    return P.BOND_T[(_ty(mol, i, j),)][0] * NM


def _t0(mol, a, b, c):  # rad
    return np.radians(P.ANGLE_T[(_ty(mol, a, b, c),)][0])


def _add(system, force, term):
    force.setForceGroup(TERM_GROUP[term])
    system.addForce(force)
    return force


def build_system(mol):
    """Return an OpenMM System holding all nine CFF93 valence terms for `mol`."""
    system = mm.System()
    for t in mol.types:
        system.addParticle(12.011 if t == "C" else 1.008)

    # --- class I part: anharmonic bond and angle, Fourier torsion -------------
    f = mm.CustomBondForce("K2*dr^2 + K3*dr^3 + K4*dr^4; dr = r - b0")
    for p in ("b0", "K2", "K3", "K4"):
        f.addPerBondParameter(p)
    for i, j in mol.bonds:
        b0, k2, k3, k4 = P.BOND_T[(_ty(mol, i, j),)]
        s = KJ / NM**2
        f.addBond(i, j, [b0 * NM, k2 * s, k3 * KJ / NM**3, k4 * KJ / NM**4])
    _add(system, f, "bond")

    f = mm.CustomAngleForce("K2*dt^2 + K3*dt^3 + K4*dt^4; dt = theta - t0")
    for p in ("t0", "K2", "K3", "K4"):
        f.addPerAngleParameter(p)
    for a, b, c in mol.angles():
        t0, k2, k3, k4 = P.ANGLE_T[(_ty(mol, a, b, c),)]
        f.addAngle(a, b, c, [np.radians(t0), k2 * KJ, k3 * KJ, k4 * KJ])
    _add(system, f, "angle")

    f = mm.CustomTorsionForce("K1*(1-cos(theta)) + K2*(1-cos(2*theta)) + K3*(1-cos(3*theta))")
    for p in ("K1", "K2", "K3"):
        f.addPerTorsionParameter(p)
    for a, b, c, d in mol.torsions():
        f.addTorsion(a, b, c, d, [k * KJ for k in P.TORSION_T[(_ty(mol, a, b, c, d),)]])
    _add(system, f, "torsion")

    # --- cross terms: all CustomCompoundBondForce ----------------------------
    # bond-bond: bonds p1-p2 and p2-p3
    f = mm.CustomCompoundBondForce(
        3, "Kbb*(distance(p1,p2)-b1)*(distance(p2,p3)-b2)"
    )
    for p in ("Kbb", "b1", "b2"):
        f.addPerBondParameter(p)
    for a, b, c in mol.bond_pairs():
        key = tuple(sorted((_ty(mol, a, b), _ty(mol, b, c))))
        f.addBond([a, b, c], [P.BOND_BOND_T[key] * KJ / NM**2, _b0(mol, a, b), _b0(mol, b, c)])
    _add(system, f, "bond-bond")

    # bond-angle: bond p1-p2 and angle p1-p2-p3 (the bond is one arm of the angle)
    f = mm.CustomCompoundBondForce(3, "Kba*(distance(p1,p2)-b0)*(angle(p1,p2,p3)-t0)")
    for p in ("Kba", "b0", "t0"):
        f.addPerBondParameter(p)
    for a, b, c in mol.angles():
        akey = _ty(mol, a, b, c)
        for i, j, k in ((a, b, c), (c, b, a)):  # arm i-j, then the other arm k
            K = P.BOND_ANGLE_T[(_ty(mol, i, j), akey)]
            f.addBond([i, j, k], [K * KJ / NM, _b0(mol, i, j), _t0(mol, a, b, c)])
    _add(system, f, "bond-angle")

    # angle-angle: angle(p1,p2,p3) and angle(p4,p2,p3); common bond p2-p3
    f = mm.CustomCompoundBondForce(
        4, "Kaa*(angle(p1,p2,p3)-t1)*(angle(p4,p2,p3)-t2)"
    )
    for p in ("Kaa", "t1", "t2"):
        f.addPerBondParameter(p)
    for a, v, common, c in mol.angle_pairs():
        k1, k2 = _ty(mol, a, v, common), _ty(mol, c, v, common)
        K = P.ANGLE_ANGLE_T[(tuple(sorted((k1, k2))), _ty(mol, v, common))]
        f.addBond([a, v, common, c], [K * KJ, _t0(mol, a, v, common), _t0(mol, c, v, common)])
    _add(system, f, "angle-angle")

    # torsion-based cross terms, particles p1..p4 = the torsion atoms a,b,c,d
    four = "K1*cos(phi)+K2*cos(2*phi)+K3*cos(3*phi)"
    f = mm.CustomCompoundBondForce(
        4,
        f"(distance(p2,p3)-bc)*({four})"
        f" + (distance(p1,p2)-bt1)*(L1*cos(phi)+L2*cos(2*phi)+L3*cos(3*phi))"
        f" + (distance(p3,p4)-bt2)*(M1*cos(phi)+M2*cos(2*phi)+M3*cos(3*phi));"
        " phi = dihedral(p1,p2,p3,p4)",
    )
    for p in ("bc", "K1", "K2", "K3", "bt1", "L1", "L2", "L3", "bt2", "M1", "M2", "M3"):
        f.addPerBondParameter(p)
    for a, b, c, d in mol.torsions():
        tkey = _ty(mol, a, b, c, d)
        K = P.BOND_TORSION_CENTRAL_T[(_ty(mol, b, c), tkey)]
        L = P.BOND_TORSION_TERMINAL_T[(_ty(mol, a, b), tkey)]
        M = P.BOND_TORSION_TERMINAL_T[(_ty(mol, c, d), tkey)]
        f.addBond(
            [a, b, c, d],
            [_b0(mol, b, c), *[k * KJ / NM for k in K],
             _b0(mol, a, b), *[k * KJ / NM for k in L],
             _b0(mol, c, d), *[k * KJ / NM for k in M]],
        )
    _add(system, f, "bond-torsion")

    f = mm.CustomCompoundBondForce(
        4,
        "(angle(p1,p2,p3)-t1)*(K1*cos(phi)+K2*cos(2*phi)+K3*cos(3*phi))"
        " + (angle(p2,p3,p4)-t2)*(L1*cos(phi)+L2*cos(2*phi)+L3*cos(3*phi));"
        " phi = dihedral(p1,p2,p3,p4)",
    )
    for p in ("t1", "K1", "K2", "K3", "t2", "L1", "L2", "L3"):
        f.addPerBondParameter(p)
    for a, b, c, d in mol.torsions():
        tkey = _ty(mol, a, b, c, d)
        K = P.ANGLE_TORSION_T[(_ty(mol, a, b, c), tkey)]
        L = P.ANGLE_TORSION_T[(_ty(mol, b, c, d), tkey)]
        f.addBond(
            [a, b, c, d],
            [_t0(mol, a, b, c), *[k * KJ for k in K], _t0(mol, b, c, d), *[k * KJ for k in L]],
        )
    _add(system, f, "angle-torsion")

    f = mm.CustomCompoundBondForce(
        4, "Kaat*(angle(p1,p2,p3)-t1)*(angle(p2,p3,p4)-t2)*cos(dihedral(p1,p2,p3,p4))"
    )
    for p in ("Kaat", "t1", "t2"):
        f.addPerBondParameter(p)
    for a, b, c, d in mol.torsions():
        k1, k2 = _ty(mol, a, b, c), _ty(mol, b, c, d)
        K = P.ANGLE_ANGLE_TORSION_T[(tuple(sorted((k1, k2))), _ty(mol, a, b, c, d))]
        f.addBond([a, b, c, d], [K * KJ, _t0(mol, a, b, c), _t0(mol, b, c, d)])
    _add(system, f, "angle-angle-torsion")
    return system


def term_energies(mol):
    """OpenMM energy of each term, in kcal/mol, evaluated on mol.xyz."""
    system = build_system(mol)
    ctx = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName("Reference"))
    ctx.setPositions(mol.xyz)
    out = {}
    for term, g in TERM_GROUP.items():
        e = ctx.getState(getEnergy=True, groups={g}).getPotentialEnergy()
        out[term] = e.value_in_unit(e.unit) / KJ
    return out


@click.command()
@click.argument("molecule", default="ethane", type=click.Choice(sorted(MOLECULES)))
def main(molecule):
    """Print each CFF93 term (OpenMM vs. plain numpy) for a perturbed MOLECULE."""
    mol = MOLECULES[molecule]()
    omm, hand = term_energies(mol), ref.energies(mol)
    click.echo(f"{molecule}: energies in kcal/mol at a randomly displaced geometry")
    click.echo(f"{'term':22s}{'OpenMM':>12s}{'numpy':>12s}{'diff':>10s}")
    for t in ref.TERMS:
        click.echo(f"{t:22s}{omm[t]:12.5f}{hand[t]:12.5f}{omm[t]-hand[t]:10.1e}")


if __name__ == "__main__":
    main()
