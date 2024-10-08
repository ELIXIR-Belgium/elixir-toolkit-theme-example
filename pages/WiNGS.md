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

## Test Connection

Save the following script as WiNGS_api.py

```python
#!/usr/bin/env python

import requests
import json
import sys
import os
import re
import pandas as pd
import os
import re
import xlsxwriter
import getpass
import pprint

# basic api interaction function
def get_api_result(url, endpoint, token, method='GET', **kwargs):
   url = url.rstrip('/')
   # add authentication token
   headers = {"Content-Type": "application/json", "Authorization": f"Bearer {token}"}
   try:
      response = requests.request(method, f"{url}/{endpoint}", headers=headers, json=arguments, params=arguments, verify=True)
      response.raise_for_status()
      return response.json()['message']
   except requests.exceptions.RequestException as e:
      print('Failed to get api response')
      print(repr(e))
      sys.exit(1)

if __name__ == '__main__':
   # provide your token : 
   token = getpass.getpass('Enter your token : ')

   url="https://wings.esat.kuleuven.be/rest-api/"
   # get the principle investigator ID (piid) associated to your user
   endpoint='piid'
   r = get_api_result(url,endpoint,token)
   print(f"The resulting data is : \n {pprint.pformat(r)}")
   # select first (if multiple)
   piid = r[0]['PIID']
```

When executed, the above script will ask for your token and will print the api response to screen. Note the piid as you will need the value later on. 

### Listing Samples

Samples access is restricted by account.  This means that as a registered user, you can only see samples assigned to you. This type of work is not federated in a strict sense, as WiNGS automatically directs queries to the correct, individual data node.  We'll move on to cross-node routines in later sections. 

For now, we'll investigate what samples you have access to, select a family (inddex + parents) and perform variant querying on them.

Add the following snippet to the __main__ section of the script:

```python
   # list all individuals
   endpoint='individuals'
   r = get_api_result(url,endpoint,token)
   print(f"You have access to {len(r)} individuals:")
   demo_sample = None
   for i in r:
      if i['LocalID'] == 'Demo_index':
         demo_sample = i
      print(f"{i['LocalID']} : {i['IndividualID']}")



```
It will list all individuals and then select "Demo_index" as an example. It is part of a public dataset all users have access to.





=============================================================================================================

The file `genesdictRefGeneAnnovarComplete_withoutzerogenes.rtf` contains a dictionary of genes, which serves as a reference for identifying whether a gene being processed is present in the list or not. An example excerpt from the file is:

```plaintext
{
  0: 'IQGAP1', 1: 'FFAR1', 2: 'ZNF829', 3: 'CNIH3', 4: 'NKD1',
  5: 'TSTD2', 6: 'DHX33', 7: 'HRK', 8: 'LOC100505841', 9: 'SPANXA2',
  10: 'GAGE12C', 11: 'TRAPPC2', 12: 'WASHC1', 13: 'APP', 14: 'CASP10',
  15: 'HECW2', 16: 'PRDM15', 17: 'SNORD25', 18: 'ATP2A2', 19: 'ZNF385A',
  20: 'PARP1', 21: 'PHF20L1', 22: 'CTIF', 23: 'ZNF561-AS1', 24: 'TMEM233',
  25: 'ITPR2', 26: 'ITGB3BP', 27: 'PROZ', 28: 'MIR8072'
}
```

This dictionary is used to verify the presence of genes during the data processing steps.



---

## Repository Contents

The repository contains the following files and directories essential for running the experiments:

- `standaloneFL.py`: Script for running Experiment 2.
- `flClient.py`: Federated Learning (FL) client script for Experiment 1.
- `flServer.py`: FL server script for Experiment 1.
- `NXTfusion`: Folder containing core source code necessary for running NXTppi.
- `marshalledP3`, `phenopediaCrohnGenesmodels`: Folders containing data needed by the scripts.
- `sources`: Additional functions required by the scripts.

---

## Running Experiment 1 (Exp1 Fedrated)

Experiment 1 involves running FedCrohn in a client-server setup with three datasets. Follow these steps:

1. **Prepare the Datasets**
   - The `marshalledP3` folder contains three CD cases control datasets (2, 3, and 4). These datasets are processed and binarized versions derived from the CAGI challenges.

2. **Open Three Bash Shells**
   - Open three terminal windows or bash shells.
   - Activate the `FedCrohn` environment in all of them by typing:
     ```sh
     conda activate FedCrohn
     ```

3. **Launch the Server**
   - In the first shell, start the FL server using dataset 4:
     ```sh
     python flServer.py 4
     ```

4. **Launch the Clients**
   - In the second shell, start the first client using dataset 3:
     ```sh
     python flClient.py 3
     ```
   - In the third shell, start the second client using dataset 2:
     ```sh
     python flClient.py 2
     ```
     
## 5. Initiating the Simulation

### Experiment 1 Workflow (Fedrated)

- The FL server (`flServer.py`) waits for the clients (`flClient.py`) to connect.
- Once connected, the server initializes an empty neural network (NN) model and sends it to the clients.
- The clients train the model on their respective datasets and send the updated models back to the server.
- The server averages these models to create a consensus model.
- This process constitutes one *round* of federated learning (FL) training. The script repeats this for five rounds and prints the final validation performance.

Feel free to experiment with different dataset combinations and observe the results.

---


##DEMO VIDEO

[![Alt text for the video](https://img.youtube.com/vi/SUPriLcMggk/0.jpg)](https://www.youtube.com/watch?v=SUPriLcMggk)