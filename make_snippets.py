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
PREAMBLE = r"""\documentclass[12pt,border=6pt,varwidth=20cm]{standalone}
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
    """Stacked layout, 12 pt text: name and OpenMM class on one line, equation below."""
    blocks = []
    for term, (hex_, family, eq, omm, ref) in S.TERMS.items():
        name = S.macro_name(term)
        if r"\left[" in eq:  # too wide for one line at 12 pt: break after the opening bracket
            head_, tail_ = eq.split(r"\left[", 1)
            tail_ = tail_.replace(r"\right]", "")
            eq = r"\begin{array}{l}" + head_ + r"[\\ \quad " + tail_ + r"\,]\end{array}"
        blocks.append(
            rf"\noindent\parbox{{16cm}}{{\textcolor{{{name}}}{{\textbf{{{term.replace('-', '--')}}}}}"
            rf"\hfill\textcolor{{{name}}}{{\texttt{{{omm}}}}}\\[1pt]"
            "\n"
            rf"\hspace*{{1em}}$\textcolor{{{name}}}{{\displaystyle {eq}}}$}}\\[7pt]"
        )
    head = (
        r"\noindent\parbox{16cm}{\textbf{Term}\hfill\textbf{OpenMM class}\\"
        r"\hspace*{1em}\textbf{Energy of one summand (Hwang et al. 1994, Eq.\ 1)}\\[-2pt]\hrule}\\[7pt]"
    )
    return (
        PREAMBLE.replace("varwidth=20cm", "varwidth=16cm")
        + head + "\n" + "\n".join(blocks) + "\n\\end{document}\n"
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
