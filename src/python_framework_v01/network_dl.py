from tqdm import tqdm
import requests

# TODO: Suppress progress bar re-print after program output
def download(fullpath, link):
    output_to_file = "{fullpath}".format(fullpath=fullpath)
    with requests.get(link, stream=True) as req:
        total_size = int(req.headers.get("content-length", 0))
        block_size = 1024
        t = tqdm(total=total_size, unit="iB", unit_scale=True, unit_divisor=block_size)
        req.raise_for_status()
        with open(output_to_file, "wb") as outfile:
            for chunk in req.iter_content(chunk_size=block_size * 4):
                t.update(len(chunk))
                if chunk:
                    outfile.write(chunk)
    return output_to_file
