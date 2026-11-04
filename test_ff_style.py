from pathlib import Path

import class2_reference as ref
import ff_style

HERE = Path(__file__).parent


def test_palette_covers_every_implemented_term():
    assert set(ff_style.TERMS) == set(ref.TERMS)


def test_colours_are_distinct():
    hexes = [v[0] for v in ff_style.TERMS.values()]
    assert len(set(hexes)) == len(hexes)


def test_chapter_colours_in_sync():
    assert (HERE / "chapter" / "colors.tex").read_text() == ff_style.latex_definitions()
