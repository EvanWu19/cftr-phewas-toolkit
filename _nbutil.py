"""
_nbutil.py — tiny helper so every toolkit notebook is built identically.

Usage:
    from _nbutil import md, code, write_notebook
    write_notebook("notebooks/01_gnomad.ipynb", [
        md("# Title\n\nSome explanation."),
        code("import toolkit as tk\ntk.load_gnomad_missense().head()"),
    ])

Notebooks are written un-executed; run them with:
    jupyter nbconvert --to notebook --execute --inplace notebooks/*.ipynb
"""
from __future__ import annotations
import nbformat as nbf


def md(text: str) -> nbf.NotebookNode:
    return nbf.v4.new_markdown_cell(text)


def code(text: str) -> nbf.NotebookNode:
    return nbf.v4.new_code_cell(text)


def write_notebook(path: str, cells: list) -> None:
    nb = nbf.v4.new_notebook()
    nb["cells"] = cells
    nb["metadata"] = {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python"},
    }
    with open(path, "w", encoding="utf-8") as fh:
        nbf.write(nb, fh)
    print("wrote", path, f"({len(cells)} cells)")
