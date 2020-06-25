import xarray as xr
import pandas as pd
import json
import os

root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def _handle_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-n",
        "--gnis-name",
        help="Enter a string corresponding to the GNIS name you wish to search for, e.g., 'Mississippi River'.",
        # TODO: accept multiple or a Path (argparse Action perhaps)
        # action='append',
        # nargs=1,
        dest="name_to_find",
        default="South Platte River",
    )

    return parser.parse_args()

def main():
    full_nc = os.path.join(
        root,
        r"test",
        r"input",
        r"geo",
        r"Channels",
        r"RouteLink_NHDPLUS.nwm.v2.0.2.nc",
    )

    names_json = os.path.join(
        root, r"test", "input", "json", r"nwm_reaches_conus_v21_wgnis_name.json",
    )

    with open(names_json, "r") as json_file:
        names = json.load(json_file)

    ids_wnames = {}
    for name, ids in names.items():
        for id in ids:
            ids_wnames[id] = name

    ids_wnames_df = pd.DataFrame.from_dict(ids_wnames, orient="index", columns=["name"])

    from networkbuilder import get_down_connections

    key_col = "link"
    downstream_col = "to"
    length_col = "Length"

    rows = (xr.open_dataset(full_nc)).to_dataframe()
    rows = rows.set_index([key_col])

    connections = get_down_connections(
        rows=rows,
        key_col=key_col,
        mask_set=set(rows.index.values),
        downstream_col=downstream_col,
        length_col=length_col,
        verbose=True,
        debuglevel=-2,
    )

    print(len(names))
    print(len(ids_wnames))
    print(len(rows))

    """ 
    Goal: 
    Find named segments tributary to the Mississippi River
    Pseudocode:

    Make a set of Named Mississippi segments
    Make a set of connections pointing to the Named Mississippi segments
    Drop from that set connections that are in the Mississippi themselves
    Also drop from that set connections that don't have a name themselves

    """

    named_segs = names["Mississippi River"]
    # named_segs = names["Missouri River"]
    print(len(named_segs))
    to_named_segs = set()
    for k, v in connections.items():
        if v[downstream_col] in named_segs:
            to_named_segs.add(k)

    tribs = to_named_segs - set(named_segs)
    named_tribs = tribs & set(ids_wnames_df.index)
    print(named_tribs)
    named_tribs_wnames = {trib: ids_wnames[trib] for trib in named_tribs}
    print(named_tribs_wnames)
if __name__ == "__main__":
    main()
