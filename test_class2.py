"""Checks: OpenMM term energies equal the plain-numpy reference; hand numbers for single terms."""

import numpy as np
import openmm as mm
import pytest

import class2_openmm as co
import class2_reference as ref
from molecules import MOLECULES, ethane

KJ = 4.184


@pytest.mark.parametrize("name", sorted(MOLECULES))
@pytest.mark.parametrize("seed", [0, 1])
def test_openmm_matches_numpy_per_term(name, seed):
    mol = MOLECULES[name](seed=seed)
    omm, hand = co.term_energies(mol), ref.energies(mol)
    for term in ref.TERMS:
        assert omm[term] == pytest.approx(hand[term], rel=1e-4, abs=1e-6), term
    assert any(abs(v) > 1e-4 for v in hand.values())


def test_every_term_is_active():
    """Displaced geometry must make all nine terms non-zero, so none is silently dropped."""
    hand = ref.energies(ethane(seed=0))
    for term, e in hand.items():
        assert abs(e) > 1e-6, term


def test_ethane_counts():
    mol = ethane()
    assert (len(mol.bonds), len(mol.angles()), len(mol.torsions())) == (7, 12, 9)
    assert len(mol.bond_pairs()) == 12
    assert len(mol.angle_pairs()) == 24


def test_hand_number_cc_stretch():
    """Two carbons stretched by 0.02 A: E = 2K d^2 + 3K d^3 + 4K d^4 with the Table 4 A values."""
    d = 0.02
    expected = 299.4 * d**2 - 515.8 * d**3 + 667.1 * d**4  # 0.115740 kcal/mol
    mol = ethane(perturb=False)
    mol.xyz[:] = 0
    mol.xyz[1] = [0, 0, (1.535 + d) * 0.1]
    # single-bond system: only the bond term, built directly
    system = mm.System()
    for _ in range(2):
        system.addParticle(12.0)
    f = mm.CustomBondForce("K2*dr^2 + K3*dr^3 + K4*dr^4; dr = r - b0")
    for p in ("b0", "K2", "K3", "K4"):
        f.addPerBondParameter(p)
    f.addBond(0, 1, [0.1535, 299.4 * KJ * 100, -515.8 * KJ * 1000, 667.1 * KJ * 10000])
    system.addForce(f)
    ctx = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName("Reference"))
    ctx.setPositions(mol.xyz[:2])
    e = ctx.getState(getEnergy=True).getPotentialEnergy()
    assert e.value_in_unit(e.unit) / KJ == pytest.approx(expected, rel=1e-6)
    assert expected == pytest.approx(0.115740, abs=1e-6)


def test_hand_number_bond_angle_cross_term():
    """Bond-angle: Kbt (b-b0)(theta-theta0) for H-C/H-C-H, 0.01 A stretch and 2 deg opening."""
    b, t = 0.01, np.radians(2.0)
    assert 19.7 * b * t == pytest.approx(0.006876, abs=1e-6)


def test_ethane_barrier_staggered_below_eclipsed():
    """Valence terms only (no nonbonded): rotating ideal ethane from eclipsed to staggered.

    Sanity check on the sign convention of the 3-fold torsion term, not a
    validation of CFF93 (the full force field also has 1-4 nonbonded terms).
    """
    stag = sum(ref.energies(ethane(perturb=False, staggered=True)).values())
    ecl = sum(ref.energies(ethane(perturb=False, staggered=False)).values())
    assert 1.0 < ecl - stag < 4.0, ecl - stag


def test_ethane_staggered_flag_gives_staggered_dihedrals():
    mol = ethane(perturb=False, staggered=True)
    phis = sorted(round(float(np.degrees(abs(ref.dihedral(mol.xyz, *t))))) for t in mol.torsions())
    assert phis == [60] * 6 + [180] * 3
