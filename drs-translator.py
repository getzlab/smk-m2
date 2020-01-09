

import requests
import subprocess


def get_gs_url(pid="STAD_BR-8592", type_="normal", url=""):
    if url[0:2] == "gs":
        gs_url =  url
    elif url[0:3] == "drs":
        token = subprocess.check_output('gcloud auth print-access-token', shell=True) 
        headers = {
            'authorization': 'bearer {}'.format(token.decode().strip()) ,
            'content-type': 'application/json'
        }
        data = '{ "url": "'+ url + '" }'
        response = requests.post('https://us-central1-broad-dsde-prod.cloudfunctions.net/martha_v2', headers=headers, data=data).json()
        #print(response)
        # [f(x) for x in sequence if condition]
        gs_url = [u_["url"] for u_ in response["dos"]["data_object"]["urls"] if u_["url"][0:2] == "gs" ][0]

    file_name = "{}/{}_url.txt".format(pid, type_)
    f1 = open(file_name, "w+")
    f1.write(gs_url)
    f1.close()

    return gs_url
    
#print(snakemake.params.normal)
get_gs_url(snakemake.wildcards.pid, "normal", snakemake.params.normal)
get_gs_url(snakemake.wildcards.pid, "tumor", snakemake.params.tumor)