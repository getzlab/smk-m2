import requests

with open(snakemake.param.credential) as json_file: 
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

normal_url = get_signed_url_from_uuid(snakemake.params.pid, 
                                      snakemake.params.normal, 
                                      id = "normal")
tumor_url= get_signed_url_from_uuid(snakemake.params.pid, 
                                    snakemake.params.tumor,
                                    id = "tumor")

