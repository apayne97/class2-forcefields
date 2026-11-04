"""One colour per term, defined once; used by the snippets, the mapping table and the LaTeX chapter.

Two families:
  * class I valence terms (what almost every force field has): cool colours
  * class II additions (anharmonic cubic/quartic parts live in the same row as
    their class I term; the cross terms are the new ones): warm colours

`latex_definitions()` prints the \\definecolor lines that chapter/colors.tex holds;
`python -m pytest test_ff_style.py` checks that the file is in sync with this module.
Colour is used in displayed equations, tables and figures only, never in running text.
"""

INK = "212745"  # text colour for everything that is not a term

# name -> (hex, family, latex equation of one summand, OpenMM class, Hwang 1994 reference)
TERMS = {
    "bond": (
        "1F5FBF", "class I (anharmonic in class II)",
        r"{}^2K_b(b-b_0)^2+{}^3K_b(b-b_0)^3+{}^4K_b(b-b_0)^4",
        "CustomBondForce", "Eq. 1, Table 4 A",
    ),
    "angle": (
        "00899B", "class I (anharmonic in class II)",
        r"{}^2K_\theta(\theta-\theta_0)^2+{}^3K_\theta(\theta-\theta_0)^3+{}^4K_\theta(\theta-\theta_0)^4",
        "CustomAngleForce", "Eq. 1, Table 4 B",
    ),
    "torsion": (
        "3949AB", "class I",
        r"{}^1K_\phi(1-\cos\phi)+{}^2K_\phi(1-\cos2\phi)+{}^3K_\phi(1-\cos3\phi)",
        "CustomTorsionForce", "Eq. 1, Table 4 C",
    ),
    "bond-bond": (
        "E4572E", "class II cross term",
        r"K_{bb'}(b-b_0)(b'-b_0')",
        "CustomCompoundBondForce (3 atoms)", "Eq. 1, Table 4 F",
    ),
    "bond-angle": (
        "E0A100", "class II cross term",
        r"K_{b\theta}(b-b_0)(\theta-\theta_0)",
        "CustomCompoundBondForce (3 atoms)", "Eq. 1, Table 4 G",
    ),
    "angle-angle": (
        "C2185B", "class II cross term",
        r"K_{\theta\theta'}(\theta-\theta_0)(\theta'-\theta_0')",
        "CustomCompoundBondForce (4 atoms)", "Eq. 1, Table 4 H",
    ),
    "bond-torsion": (
        "8C5A2B", "class II cross term",
        r"(b-b_0)\left[{}^1K_{\phi b}\cos\phi+{}^2K_{\phi b}\cos2\phi+{}^3K_{\phi b}\cos3\phi\right]",
        "CustomCompoundBondForce (4 atoms)", "Eq. 1, Table 4 J, K",
    ),
    "angle-torsion": (
        "3E8E41", "class II cross term",
        r"(\theta-\theta_0)\left[{}^1K_{\phi\theta}\cos\phi+{}^2K_{\phi\theta}\cos2\phi+{}^3K_{\phi\theta}\cos3\phi\right]",
        "CustomCompoundBondForce (4 atoms)", "Eq. 1, Table 4 L",
    ),
    "angle-angle-torsion": (
        "8E24AA", "class II cross term",
        r"K_{\phi\theta\theta'}(\theta-\theta_0)(\theta'-\theta_0')\cos\phi",
        "CustomCompoundBondForce (4 atoms)", "Eq. 1, Table 4 I",
    ),
}


def macro_name(term):
    return "c" + "".join(w.capitalize() for w in term.split("-"))


def latex_definitions():
    lines = [f"\\definecolor{{ffInk}}{{HTML}}{{{INK}}}"]
    for term, (hex_, family, *_rest) in TERMS.items():
        lines.append(f"\\definecolor{{{macro_name(term)}}}{{HTML}}{{{hex_}}}  % {term}: {family}")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    print(latex_definitions(), end="")
