import requests
import json

import argparse


# run once for normal, another time for matched tumor
parser = argparse.ArgumentParser(description= 'localize {json} {uuid} {pid} {type}')

parser.add_argument('credential_json', action="store")
parser.add_argument('uuid', action="store")
parser.add_argument('pid', action="store")
parser.add_argument('type', action="store")
result = parser.parse_args()

# get the uuid for indices
def get_bai_uuid(uuid="736a8e90-85ec-4007-b34a-1bf823eec6fc", id="normal"):
    file_endpt = 'https://api.gdc.cancer.gov/legacy/files/'
    file_with_indice = '?expand=index_files'
    response = requests.get(file_endpt + uuid + file_with_indice)
    res = response.json()
    #return(res)
    bai_uuid = res['data']['index_files'][0]['file_id']
    return(bai_uuid)

get_bai_uuid()

# start with NCI data commons
with open("../"+result.credential_json) as json_file: 
         credential = json.load(json_file)

token = requests.post('https://nci-crdc.datacommons.io/user/credentials/api/access_token', 
        json=credential).json()
headers = {'Authorization': 'bearer '+ token['access_token']}



def get_signed_url_from_uuid(pid, uuid = "736a8e90-85ec-4007-b34a-1bf823eec6fc",id = "normal"):
    file_endpt = 'https://nci-crdc.datacommons.io/user/data/download/'
    response = requests.get(file_endpt + uuid, headers=headers)
    url = response.json()['url']
    ### write the url to a txt file
    filename=pid+"/" + id + "_url.txt"
    file1 = open(filename,"w")
    file1.write(url)
    file1.close() 
    return(url)


index_uuid = get_bai_uuid(result.uuid, result.type)
url = get_signed_url_from_uuid(result.pid, result.uuid, result.type)
url = get_signed_url_from_uuid(result.pid, index_uuid, result.type + "_index")
