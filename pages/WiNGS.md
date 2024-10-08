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
      return response.json()
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


```

When executed, the above script will ask for your token and will print the api response to screen. Note the piid as you will need the value later on. 




## Downloading and Installing ANNOVAR

1. **Download & Install ANNOVAR**: Download the appropriate Annovar version for your operating system from the official website [here](https://www.openbioinformatics.org/annovar/annovar_download_form.php). Follow the installation instructions provided.
2. **Download & Install reference databases**: Annovar utilizes various reference databases for annotations. Download the relevant databases based on your research focus (e.g., gene models, functional annotations). Refer to the [Annovar documentation](https://annovar.openbioinformatics.org/en/latest/user-guide/startup/) for downloading and installing these databases in the designated directory.
3. **Update Databases**: Download the necessary annotation databases to ensure ANNOVAR has the latest information for your analysis.

### Download the Databases

1. Visit the [Annovar download page](https://annovar.openbioinformatics.org/en/latest/user-guide/download/). This page lists available databases for download, including human reference databases.
2. **Download refGene database**: This database provides information about protein-coding genes. Click on the "download" button next to the refGene database. This will typically download a compressed archive file (e.g., .gz).
3. **Create the humandb Directory**: Navigate to the directory where you want to store your Annovar databases. Create a new folder named `humandb` using the appropriate command for your operating system (e.g., `mkdir humandb` in Linux/macOS).
4. **Extract Downloaded Archives**: Move the downloaded archive files for each database into the `humandb` directory. Use an archive extraction tool like `gunzip` (Linux/macOS) or a dedicated archive manager to extract the contents of each file.

### Configuring ANNOVAR

1. **Set Input File**: Specify the VCF file you want to annotate using ANNOVAR.
2. **Choose Databases**: Select the appropriate annotation databases to use for your analysis (e.g., dbNSFP version 42a).
3. **Customize Parameters**: Adjust ANNOVAR settings to suit your specific research needs and preferences.

### Annotating the VCF File with ANNOVAR

1. **Run ANNOVAR**: Execute the ANNOVAR command to start the annotation process using the `table_annovar.pl` script provided with Annovar. Here's the basic command structure:
   ```sh
   perl .../annovar/table_annovar.pl VCFFILE .../annovar/humandb/ -out OUTPUT -vcfinput -buildver hg19 -protocol refGene,dbnsfp42a -operation gxf -xreffile XREFFILE
   
   
## 3. Converting Genetic Annotations to Vector Representations

This section details the data processing steps performed for the FedCrohn project. The process involves working with gene dictionaries, label processing, and feature extraction from annotated VCF files.

#### 1. Gene Dictionary

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

#### 2. Label Processing

The labels for the samples are processed based on their identifiers. The code snippet for processing the labels is as follows:

```plaintext
if row['actual_sample'].startswith('CD'):
    labels.loc[idx, 'actual_sample'] = 0
else:
    labels.loc[idx, 'actual_sample'] = 1
```

In this context:
- Samples with identifiers starting with 'CD' are labeled as `0`.
- All other samples are labeled as `1`.

This labeling helps distinguish between different types of samples in the dataset.

#### 3. Feature Extraction from Annotated VCF Files

The core of the data processing involves iterating over each of the processed files, passing them to the `parseAnnovarIntoFeatures` method. This method reads Variant Call Format (VCF) files annotated by ANNOVAR and creates a matrix containing information on variants per gene. The matrix is returned as a NumPy ndarray.

The steps involved are as follows:

1. **Read Annotated VCF Files**: The VCF files, annotated by ANNOVAR, contain detailed information about genomic variants.
2. **Create a Matrix of Variants**: The `parseAnnovarIntoFeatures` function parses these files and constructs a matrix where each entry corresponds to a variant in a specific gene.
3. **Store Feature Vectors**: The resulting matrices (feature vectors) are stored in a list of genomes.
4. **Save Feature Vectors**: These feature vectors are saved as `.npy` files, which serve as input for subsequent analysis and modeling.

The `parseAnnovarIntoFeatures` method performs several key tasks:
- **File Reading and Validation**: Opens the VCF file, reads it line by line, and validates the presence of required columns (e.g., "Chr", "Start").
- **Region and Gene Filtering**: Filters out irrelevant regions and genes not present in the reference gene list.
- **Variant Type Handling**: Handles various types of variants, including synonymous and nonsynonymous SNVs, and assigns them to the correct positions in the matrix.
- **Frequency Thresholding**: Applies frequency thresholds to variants if a frequency dictionary is provided.
- **Scoring**: Assigns functional scores to variants based on different predictors (e.g., provean, metasvm) if applicable.

#### Summary of Data Processing Steps

1. **Gene Dictionary Setup**: Use the provided gene dictionary to validate genes during processing.
2. **Label Processing**: Assign numerical labels to samples based on their identifiers.
3. **Feature Extraction**: Parse annotated VCF files to create matrices of gene variants, filtering and scoring variants as necessary.
4. **Saving Features**: Save the processed feature vectors as `.npy` files for further analysis.

By following these steps, we ensure that the data is processed consistently and is ready for use in the federated learning setup of FedCrohn. This involves leveraging a gene dictionary for validation, processing sample labels, and extracting meaningful features from annotated genomic data.


---

## 4. Setting Up the Federated Learning Repository

To run FedCrohn, a federated learning system designed for Crohn's disease prediction, it's recommended to use Miniconda for managing the Python environment and dependencies. Follow the steps below to set up the environment:

### Prerequisites

- **Miniconda**: Miniconda is a minimal installer for Conda, which is an open-source package management and environment management system. Miniconda allows you to create isolated environments for different projects.

### Step-by-Step Setup Instructions

1. **Download and Install Miniconda**
   - Visit the [Miniconda installation page](https://docs.conda.io/en/latest/miniconda.html).
   - Download the appropriate installer for your operating system (Windows, macOS, or Linux).
   - Follow the installation instructions provided on the website.

2. **Create a New Conda Environment**
   - Open your terminal or command prompt.
   - Create a new environment named `FedCrohn` with Python 3.7 by running:
     ```sh
     conda create -n FedCrohn python=3.7
     ```

3. **Activate the Environment**
   - Activate the newly created environment with the command:
     ```sh
     conda activate FedCrohn
     ```

4. **Install PyTorch**
   - Install PyTorch version 1.0 or higher. You can use the following command or refer to the [PyTorch website](https://pytorch.org) for more options:
     ```sh
     conda install pytorch -c pytorch
     ```

5. **Install SciPy and NumPy**
   - Install SciPy and NumPy libraries with the command:
     ```sh
     conda install numpy scipy
     ```

6. **Install Flower (FL Library)**
   - Flower is the Federated Learning library required for running FedCrohn. Refer to the [Flower documentation](https://flower.dev/docs/) for installation instructions.

### Additional Notes

- During the setup process, you may encounter dependency resolutions and patches that need to be applied. Ensure you follow any prompts or error messages to install the necessary versions of libraries and resolve conflicts.
- Some of the libraries may need to be built from source rather than installed via pip. This may involve cloning the library's repository and following the build instructions provided.

### Removing the Environment

If you need to remove the `FedCrohn` environment at any time, you can do so by typing:
```sh
conda remove -n FedCrohn --all
```

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



### Running Experiment 2 (Exp2 Standalone)

Experiment 2 simulates the server-client interaction using a single script `standaloneFL.py`. This script combines all three datasets and simulates FL training with multiple clients.

1. **Run the Script**
   - Execute the following command to run the script:
     ```sh
     python standaloneFL.py 4
     ```
   - The script will automatically split the data, assign it to three centers, and simulate the FL process.

### Experiment 2 Workflow

- The script simulates the performance of FL methods for Crohn's disease prediction when several small datasets are combined into a single federated learning effort.

---

By following these instructions, you will set up the necessary environment and run the experiments for FedCrohn. If you have any questions or need further assistance, please follow up as needed.

##DEMO VIDEO

[![Alt text for the video](https://img.youtube.com/vi/SUPriLcMggk/0.jpg)](https://www.youtube.com/watch?v=SUPriLcMggk)