---
title: Tutorial on using the WiNGS REST API to query Federated Genomics Data 
---


### Base Repo :  [wings-public](https://github.com/wings-public)
### Homepage : [wings-platform.org](https://wings-platform.org)
## Introduction

This tutorial will guide you through the process of querying federated genomics data hosted in WiNGS, using python based access to the REST api. We'll cover the following sections:

1. Understanding Genetic Data and VCF Files
2. Setting up the python conda environment
3. Acquiring an api key
4. Querying data 
5. Manipulating data

## 1. Understanding Genetic Data, VCF Files and the federated setup of WiNGS

### VCF files

Genetic data is usually provided in Variant Call Format (VCF) files. These files contain information about genetic variants, including the position, reference, and alternate alleles, as well as additional annotations. VCF files are commonly used in genomic research to store and share genetic variant data.  It is important to keep in mind that over time, multiple versions of '[the human genome](https://en.wikipedia.org/wiki/Reference_genome)' have been released, and variant positions might shift slightly between these releases.  

VCF files can contain data from a single sample, or from multiple samples. In the latter case, samples are organized is columns, holding information for each encountered position. Although these files can be filtered and analyzed using multiple tools such as [vcftools](https://github.com/vcftools/vcftools) or [pyVCF](https://pyvcf.readthedocs.io/en/latest/), adding more data over time would require to rebuild them iteratively.  Another limitation is formed by the privacy constraints on genomic data, prohibiting the collection of cross-center cohorts for large-scale statistical analysis. 

### Principle of WiNGS 

For these reasons, WiNGS was implemented to become a portal to aggregated genomics data, kept privately in a federated framework. WiNGS supports data in multiple of the human reference versions (or genome builds) and provides unified, biological annotations. This way, variants in all federated data hubs, will have uniform annotations (or meta data). 

TODO: add wings setup img : how ? 

Data is structured in WiNGS as follows:
- ** Individual **: these are cinical cases, patients.
- ** Sample **: these are the experimental samples associated to individuals. (eg : WGS experiment, or SV experiment)
- ** Dataset **: These are the analysis files related to a sample. eg : GATK based SNV analysis on hg38 for WGS data. 

WiNGS is accessible as REST api, which is documented through [swagger](https://wings.esat.kuleuven.be/rest-api/api-docs/).  To perform an analysis, the following prerequisites are needed, which are performed through a web-based UI: 

- ** A user account ** : to shield data, anonymous access is not allowed, register [here](https://wings.esat.kuleuven.be/Account/Register)
- ** An API key ** : [Log in](https://wings-platform.org) and select user meny : Create API key , on the top right



 

## 2. Setup conda environment

### Install miniconda

We use miniconda for this tutorial. Follow the [Installation Instructions](https://docs.anaconda.com/miniconda/miniconda-install/). For example on linux: 

```bash
   mkdir -p ~/miniconda3
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
   bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
   rm ~/miniconda3/miniconda.sh
```

Install the prerequisites : 

```bash
   conda create -n wings  python=3.10 pandas xlsxwriter 
   conda activate wings
```

## 3. Test Connection

Save the following script as WiNGS_api.py

```python
#!/usr/bin/env python

@title IMPORTS
import getpass
import requests
import time
from tqdm import tqdm
import pprint
import pandas as pd
import matplotlib.pyplot as plt
import hashlib
import json
%matplotlib inline

@title GENERAL CLASSES
#############
## CLASSES ##
#############
class APIError(ValueError):
    def __init__(self, response, value: dict):
        self.response = response
        self.value : dict = value

class WingsApi:

    SLEEP_SECONDS  = 5
    RETRY = 60 # number of times to retry, wait SLEEP seconds after each failure

    def __init__(self,url,token):
        self.url = url.rstrip('/')
        self.token = token
        self._job_cache = {}
        self.session = requests.session()
        self.session.headers.update({'Authorization': f'Bearer {self.token}'})
        self.session.headers.update({'Content-Type': 'application/json'})

    def hash_request(self,data):
        # Convert dictionary to a JSON string with sorted keys to ensure consistent ordering
        json_str = json.dumps(data, sort_keys=True)
        # Create an MD5 hash object
        md5_hash = hashlib.md5()
        # hash the json
        md5_hash.update(json_str.encode('utf-8'))
        # Get the hexadecimal representation of the hash
        return md5_hash.hexdigest()

    def get(self, endpoint, arguments={},skip_cache=False):
        md5 = self.hash_request([endpoint, arguments])
        if md5 in self._job_cache and not skip_cache:
            return self._job_cache[md5]
        else:
          r = self.session.get(f"{self.url}/{endpoint}",params=arguments)
          if r.status_code != 200:
            raise APIError(r,r.json())
          self._job_cache[md5] = r.json()['message']
          return r.json()['message']

    def post(self,endpoint,arguments={},skip_cache=False):
        md5 = self.hash_request([endpoint, arguments])
        if md5 in self._job_cache and not skip_cache:
            return self._job_cache[md5]
        else:
          r = self.session.post(f"{self.url}/{endpoint}",json=arguments)
          if r.status_code != 200:
            raise APIError(r,r.json())
          self._job_cache[md5] = r.json()['message']
          return r.json()['message']

    def individual_details(self,individual):
         # clinical data
         endpoint=f"individuals/{individual['IndividualID']}"
         r = self.get(endpoint)
         individual['phenotype'] = r[0]['phenotype']
         # trio's defined ?
         endpoint='families'
         r = self.get(endpoint,{'proband' : individual['IndividualID']})
         if len(r):
            individual['family'] = r
            endpoint = 'trios'
            r = self.get(endpoint,{'family' : individual['family']['FamilyID']})
            r = [x for x in r if x['TrioStatus'] != 'disabled']
            if len(r):
              individual['trio'] = r
            else:
              individual['trio'] = None
         else:
            individual['family'] = None
            individual['trio'] = None
         # samples (data)
         endpoint='samples'
         r = self.get(endpoint,{'host_id' : 1, 'piid' : individual['PIID']})
         r = [x for x in r if x['IndividualID'] == individual['IndividualID']]
         if len(r):
           individual['samples'] = r
         else:
           individual['samples'] = r
         return individual

    def get_results(self,endpoint,arguments):
        results = list()
        # first call : with cache
        r = self.get(endpoint,arguments)
        if r.get('status') == 'ready':
          results.extend(r['results'])
          # pages ?
          page = 1
          while r.get('meta').get('last_page') == 'false':
            page += 1
            #print(f"Fetching page {page}...")
            arguments['page'] = page
            r = self.get(endpoint,arguments)
            results.extend(r['results'])
          return results

        # not ready / cached: poll status
        for _ in tqdm(range(self.RETRY)):
          r = api.get(endpoint,arguments,skip_cache=True)
          status = r.get("status")
          if status == "inprogress":
              time.sleep(self.SLEEP_SECONDS)
              continue
          elif status == "ready":
              page = 1
              results.extend(r['results'])
              while r.get('meta').get('last_page') == 'false':
                page += 1
                print(f"Fetching page {page}...")
                arguments['page'] = page
                r = self.get(endpoint,arguments)
                results.extend(r['results'])
          raise ValueError(r)

        print(f"\nJOB: {request_id} still in progress. Stop polling... ")
        return None



if __name__ == '__main__':
   # provide your token : 
   token = getpass.getpass('Enter your token : ')
   # the url of the api :
   url="https://wings.esat.kuleuven.be/rest-api/"
   # setup the api class
   api = WingsApi(url,token)

   # get the principle investigator ID (piid) associated to your user
   endpoint='piid'
   r = get_api_result(url,endpoint,token)
   print(f"The resulting data is : \n {pprint.pformat(r)}")
   # select first (if multiple)
   piid = r[0]['PIID']
```

When executed, the above script will ask for your token and will print the api response to screen. Note the piid as you will need the value later on. 

## 4. Listing Samples

Samples access is restricted by account.  This means that as a registered user, you can only see samples assigned to you. This type of work is not federated in a strict sense, as WiNGS automatically directs queries to the correct, individual data node.  We'll move on to cross-node routines in later sections. 

For now, we'll investigate what samples you have access to, select a family (inddex + parents) and perform variant querying on them.

Add the following snippet to the __main__ section of the script:

```python
   # list all individuals
   endpoint='individuals'
   r = api.get(endpoint)
   print(f"You have access to {len(r)} individuals:")
   individual = None
   
   # input for matching the individual you want :
   #s_string = input('Provide the LocalID of the individual to select: ')
   s_string = "Demo_index"
   
   for idx, i in enumerate(r):
      if i['LocalID'] == s_string:
         individual = i
      if idx < 5:
          print(f"{i['LocalID']} : {i['IndividualID']}")
      if idx == 5:
          print('....')
   
   if not individual:
     raise Exception('No Individual found')
   
    
   print("\nDetails on the selected individual : ")
   pprint.pprint(individual)
   
   
```
It will list all individuals and then select "Demo_index" as an example. It is part of a public dataset all users have access to.


### Information about the individual of interest

The following code will fetch information about the family members linked to the sample and set phenotypes.

```python
   # get more details, see class for details on endpoints used
   individual = api.individual_details(individual)

   print("\nDetails on the selected individual : ")
   pprint.pprint(individual)
   
   # get variables for use in querying :
   used_trio = individual['trio'][0]
   used_sample = [x for x in individual['samples'] if x['SampleFileID'] == used_trio['ProbandFileID']][0]


```

## 5. Querying Data : Sample Based Analysis

### General querying


### Family based querying


### Create resulting XLX file



## 6. Querying DAta : Variant/Gene Based Analysis

### Variant Based


### Gene Based


### 7. Querying Data : Structureal Variant 

=============================================================================================================



---


##DEMO VIDEO

[![Alt text for the video](https://img.youtube.com/vi/SUPriLcMggk/0.jpg)](https://www.youtube.com/watch?v=SUPriLcMggk)