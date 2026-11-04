# class2-forcefields

Learning repo on **class II force fields**: what the functional form is, why it has cross terms, and how each term is written in OpenMM. Everything is based on Hwang, Stockfisch & Hagler, *J. Am. Chem. Soc.* **116**, 2515 (1994) (CFF93 for alkanes); equation and table numbers in the code refer to it. Other papers in `chapter/refs.bib` are cited only for context. Not from that paper, and so not covered: PCFF, COMPASS, MM3 details.

## Setup

```bash
cd ~/science/class2-forcefields
conda env create -f environment.yaml   # creates the `class2-forcefields` env (OpenMM, numpy, click, pytest)
conda activate class2-forcefields
```

## Colours

One colour per term, defined once in `ff_style.py`. Cool colours (blue, teal, indigo) are the terms a class I force field also has; warm colours are the class II cross terms. Used in the equation snippets, the mapping table and the chapter. `python ff_style.py > chapter/colors.tex` regenerates the LaTeX colour definitions; `test_ff_style.py` checks they are in sync.

## Files

| File | What it is |
|---|---|
| `cff93.py` | CFF93 hydrocarbon parameters, transcribed from Hwang 1994 Table 4 (A–C, F–L). Nonbonded parts (D, E) not included |
| `molecules.py` | Ethane and trans n-butane, with enumerated bonds, angles, torsions and pairs; randomly displaced so every term is non-zero (not optimised) |
| `class2_reference.py` | Plain-numpy energy of each of the nine valence terms, straight from Eq. 1 |
| `class2_openmm.py` | OpenMM system with each term in its own force group (`CustomBondForce`, `CustomAngleForce`, `CustomTorsionForce`, `CustomCompoundBondForce`); `python class2_openmm.py ethane` prints OpenMM vs numpy per term |
| `test_class2.py` | OpenMM equals numpy per term (ethane, butane, two displacements); hand numbers; term counts for ethane; staggered ethane below eclipsed (valence only) |
| `ff_style.py`, `make_snippets.py` | Palette + equation text; builds `figures/term-*.{svg,png,pdf}` and `figures/openmm-mapping-table.*` |
| `chapter/` | LaTeX summary (`class2-forcefields.tex`, `refs.bib`, `colors.tex`) |

```bash
python -m pytest -q
python class2_openmm.py butane
python make_snippets.py                     # regenerate figures (needs pdflatex, pdftocairo)
cd chapter && latexmk -pdf class2-forcefields.tex
```

## PowerPoint

`figures/*.svg` (vector, Insert > Pictures) or `figures/*.png` (transparent, 300 dpi). One file per term plus the whole mapping table.

## Caveats

* Parameters are the published CFF93 alkane values, but only the valence terms are implemented, so energies are not full CFF93 energies.
* The tests show OpenMM matches an independent numpy transcription of the same formulas. Which atoms enter each sum (see the chapter) is a reading of Eq. 1 and Table 4.
* Molecules are random displacements of rough geometries, to exercise every term, not minimised structures.
