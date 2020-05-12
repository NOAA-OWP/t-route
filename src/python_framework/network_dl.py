import requests

def download(fullpath, link):
    output_to_file = "{fullpath}".format(fullpath=fullpath)
    with requests.get(link, stream=True) as req:
        req.raise_for_status()
        with open(output_to_file , "wb") as outfile:
            for chunk in req.iter_content(chunk_size=8192):
                if chunk:
                    outfile.write(chunk)
    return output_to_file

