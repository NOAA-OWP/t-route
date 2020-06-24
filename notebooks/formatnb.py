import nbformat
import black


def clean_outputs(nb):
    for cell in nb.cells:
        if cell["cell_type"] == "code":
            cell["execution_count"] = None
            cell["outputs"] = []
    return nb


def black_code(nb):
    fm = black.FileMode()
    for cellnum, cell in enumerate(nb.cells):
        if cell["cell_type"] == "code":
            try:
                cell["source"] = black.format_str(cell["source"], mode=fm)
            except black.InvalidInput:
                print(f"Possible invalid python encountered in cell {cellnum}")
                print(f"Cell source causing exception:\n{cell['source']}\n")
                continue
    return nb


def clean_metadata(nb):
    nb.metadata.language_info.version = ""
    return nb


def write_notebook(nb, filename):
    nbformat.validate(nbdict=nb)
    nbformat.write(nb, filename)


def main(args):
    with open(args.notebook, "r") as fin:
        nb = nbformat.read(args.notebook, nbformat.NO_CONVERT)

    nb = black_code(nb)
    nb = clean_outputs(nb)
    nb = clean_metadata(nb)
    write_notebook(nb, args.notebook)


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("notebook")

    args = p.parse_args()
    main(args)
