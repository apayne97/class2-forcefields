"""Build PowerPoint-ready colour-coded equation and table snippets.

For each term and for the OpenMM mapping table this writes a standalone LaTeX
file in snippets/, compiles it with pdflatex, and exports figures/<name>.pdf,
.svg (vector, via pdftocairo) and .png (transparent, 300 dpi, via pdftocairo).
Insert the SVG in PowerPoint (Insert > Pictures) or paste the PNG.
"""

import subprocess
from pathlib import Path

import click

import ff_style as S

HERE = Path(__file__).resolve().parent
PREAMBLE = r"""\documentclass[border=6pt,varwidth=20cm]{standalone}
\usepackage{amsmath,xcolor,array}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\input{../chapter/colors.tex}
\begin{document}
\color{ffInk}
"""


def equation_snippet(term):
    hex_, family, eq, omm, ref = S.TERMS[term]
    name = S.macro_name(term)
    label = term.replace("-", "--")
    return (
        PREAMBLE
        + rf"\textcolor{{{name}}}{{\textbf{{{label}}}}} \quad {{\small\color{{{name}}}{family}}}\\[2pt]"
        + "\n"
        + rf"$\displaystyle \textcolor{{{name}}}{{{eq}}}$"
        + "\n\\end{document}\n"
    )


def mapping_table_snippet():
    rows = []
    for term, (hex_, family, eq, omm, ref) in S.TERMS.items():
        name = S.macro_name(term)
        rows.append(
            rf"\textcolor{{{name}}}{{\textbf{{{term.replace('-', '--')}}}}} & "
            rf"$\textcolor{{{name}}}{{\displaystyle {eq}}}$ & "
            rf"\textcolor{{{name}}}{{\texttt{{{omm.replace('_', ' ')}}}}} \\[4pt]"
        )
    body = "\n".join(rows)
    return (
        PREAMBLE.replace("varwidth=20cm", "varwidth=30cm")
        + r"\renewcommand{\arraystretch}{1.1}"
        + "\n"
        + r"\begin{tabular}{@{}l l l@{}}"
        + "\n\\textbf{Term} & \\textbf{Energy of one summand (Hwang et al. 1994, Eq.\\ 1)} & \\textbf{OpenMM} \\\\\\hline\n"
        + body
        + "\n\\end{tabular}\n\\end{document}\n"
    )


def build(name, tex):
    (HERE / "snippets").mkdir(exist_ok=True)
    (HERE / "figures").mkdir(exist_ok=True)
    path = HERE / "snippets" / f"{name}.tex"
    path.write_text(tex)
    out = HERE / "figures"
    subprocess.run(
        ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", f"-output-directory={out}", path.name],
        cwd=path.parent, check=True, capture_output=True,
    )
    pdf = out / f"{name}.pdf"
    subprocess.run(["pdftocairo", "-svg", str(pdf), str(out / f"{name}.svg")], check=True)
    subprocess.run(
        ["pdftocairo", "-png", "-transp", "-r", "300", "-singlefile", str(pdf), str(out / name)], check=True
    )
    for ext in ("aux", "log"):
        (out / f"{name}.{ext}").unlink(missing_ok=True)


@click.command()
def main():
    """Regenerate every snippet in figures/ (SVG, PNG, PDF)."""
    for term in S.TERMS:
        build(f"term-{term}", equation_snippet(term))
        click.echo(f"term-{term}")
    build("openmm-mapping-table", mapping_table_snippet())
    click.echo("openmm-mapping-table")


if __name__ == "__main__":
    main()
