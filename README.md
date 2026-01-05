# ERASOR - Erasmus ASO designer

## Table of contents
- [🧬 Introduction](#introduction)
- [🎯 Quick start guide](#quick-start-guide)
- [📋 Features](#features)
  - [⚙️ Settings & Filters](#️-settings--filters)
    - [Gene Selection](#gene-selection)
    - [Analysis Features](#analysis-features)
    - [System Requirements](#system-requirements)
    - [Sequence Parameters](#sequence-parameters)
    - [Quality Control Thresholds](#quality-control-thresholds)
    - [ASO Modifications](#aso-modifications)
    - [Off-Target Analysis](#off-target-analysis)
  - [Page overview](#page-overview)
    - [Main page](#main-page)
    - [RNaseH page](#rnaseh-page)
    - [Off-targets page](#off-targets-page)
- [Setup guide](#setup-guide)
  - [👥 Regular user setup (Shiny server only)](#-regular-user-setup-shiny-server-only)
  - [💻 Developer setup (RStudio + Git + SSH)](#-developer-setup-rstudio--git--ssh)
  - [✔ Setup Summary Table](#-setup-summary-table)
- [📦 Dependencies](#-dependencies)


## 🧬 Introduction
Erasmus ASO designer (ERASOR) is an Antisense oligo sequence design tool for RNase H1 mediated mRNA degradation. The effectiveness of an ASO to silence a given gene is predicated on many different variables. ERASOR gives insight into what ASOs will and will not be worth testing in the lab. In silico prediction is much easier and faster than lab work. With this in mind, ERASOR is meant to be used as a first step in an ASO drug development pipeline.


## 🎯 Quick start guide
To get started quickly and get results after starting the shiny server follow these steps:

1. To start enter the **ensemble ID** for the gene you want to target. 
2. Tick or untick the checkboxes to change the output.
    - Polymorphism analysis
    - Conserved & Orthology
    - Running on Linux-OS
3. Select your desired **Oligo length** range
4. Enter your values for **filtering**
5. Press **"Run"** to start the application

The runtime of ERASOR may be 10-60 minutes. Runtime is heavily dependent on the selected oligo length and gene size and can take longer than 60 minutes if unfavorable settings are selected. 

## 📋 Features

### ⚙️ Settings & Filters

#### Gene Selection

**Ensembl Gene ID**  
Enter a valid Ensembl Gene ID (`ENSG…`) corresponding to the gene of interest.

---

#### Analysis Features

**Polymorphism Analysis**  
Identifies single nucleotide polymorphisms (SNPs) associated with the selected Ensembl gene. SNPs are genetic variations where one nucleotide is substituted for another (e.g., C → T).

**Conservation & Orthology**  
Compares the selected Ensembl gene with the mouse genome to identify conserved regions and orthologous genes. Due to the strong genetic similarity between humans and mice, this analysis helps assess whether the gene can be studied *in vivo* using mouse models.

---

#### System Requirements

**Running on Linux OS**  
Enable this option only when running the program on a Linux operating system. Certain features (e.g. ViennaRNA analysis) require Linux and may fail on Windows or macOS if enabled.

---

#### Sequence Parameters

**Oligo Length**  
Specifies the length or range of target mRNA sequences and complementary ASO sequences.  
- Larger lengths/ranges → shorter runtime  
- Smaller lengths/ranges → longer runtime

---

#### Quality Control Thresholds

**Amount of Perfect Matches**  
Defines the number of perfect matches allowed. This value is used as a quality control threshold for target mRNA sequences.

**Amount of Mismatches**  
Defines the number of mismatches allowed. This value is used as a quality control threshold for target mRNA sequences.

**Accessibility Score**  
Estimates how accessible the target mRNA sequences are, serving as a quality control metric. Higher accessibility is an important factor in identifying potentially effective ASOs.

**Polymorphism Frequency**  
Estimates the probability that a polymorphism (SNP) is present in the target mRNA sequence.  
Lower values are preferred, as they indicate a reduced likelihood of polymorphic variation.

**Toxicity Score**  
Estimates the potential toxicity of the complementary ASO targeting the mRNA sequence.  
Higher scores indicate lower toxicity and are therefore preferable.

**Motif Correlation Score**  
Certain sequence motifs are correlated with ASO effectiveness. Each motif is assigned a weight based on the strength of its correlation, and the total score is calculated by summing these weights.

---

#### ASO Modifications

**5′ Modification (5′Mod)**  
Adds chemical modifications to the 5′ end of the ASO sequence, affecting accessible binding sites on the target mRNA. An overlap of up to **2 modified nucleotides** on the 3′ end of the target mRNA is allowed.

**3′ Modification (3′Mod)**  
Adds chemical modifications to the 3′ end of the ASO sequence, affecting accessible binding sites on the target mRNA. An overlap of up to **4 modified nucleotides** on the 5′ end of the target mRNA is allowed.

---

#### Off-Target Analysis

**Off-Target Mismatches**  
Identifies off-target sequences based on the number of mismatches allowed in the target mRNA.  
  - Default: **1 mismatch**  
  - Allowing more mismatches increases off-target detection but also increases processing time.

***
### ️Page overview

#### Main page:  

This is the main page of the application, which displays the primary output table. The table lists all potential target mRNA sequences for the provided Ensembl ID and includes information to help the user select suitable ASO targets. After selecting an ASO, the application automatically redirects to the RNase H tab, and the chosen ASO becomes available for further analysis on all other tabs. 

#### RNaseH page:  

This page predicts the optimal binding site for RNase H on the chosen target mRNA sequence. It uses a matrix of dinucleotide values and a sliding window approach to evaluate all possible binding sites along the mRNA sequence. The script accounts for chemical modifications on the ASO’s 5’ and 3’ ends, which can affect possible binding sites. Each potential site is scored based on the average from the dinucleotide matrix, and the results are ranked to highlight the most promising locations for RNase H cleavage. 

#### Off-targets page  

This page displays the off-targets for the selected target mRNA sequence based on the allowed number of mismatches. For the complementary ASO sequence, it shows mismatches, deletions, insertions, protein name and additional information from GGGenome. The protein name is then given to Protein Atlas to retrieve the tissue expression data. This page also calculates the accessibility potential of each off-target.   

## Setup guide
### 👥 Regular user setup (Shiny server only)

This section is for end-users who only need to run the Shiny application. No RStudio or Git setup is required. 

#### 1. Clone git repo
Clone this git repository and navigate to the application folder.

```
git clone https://github.com/JimReijnhoudt/ASOstool-v2.git
cd ASOstool-v2
```

#### 2. Build the Docker image 

```
docker build -t asostoolv2 .
```

#### 3. Start the container


On Linux/Mac:
```         
docker run -d \
    -p 3838:3838 \
    -v $(pwd):/srv/shiny-server/ASOstool-v2 \
    asostoolv2
```

On Windows (Powershell):
```
docker run -d `
  -p 3838:3838 `
  -v ${PWD}:/srv/shiny-server/ASOstool-v2 `
  asostoolv2
```

#### 4. Open the Shiny App

Open:

👉 <http://localhost:3838/ASOstool-v2/Simpel/> 

You're done!

#### 5. Stop the container

Stop the container after you're done to free up ports.
List running containers:

```
docker ps
```

Stop it:

```
docker stop <container_id>
```

Remove it (optional):

```
docker rm <container_id>
```
***
### 💻 Developer setup (RStudio + Git + SSH)

This section is for developers contributing to the project.

#### 1. Clone git repo
Clone this git repository and navigate to the application folder.

```
git clone https://github.com/JimReijnhoudt/ASOstool-v2.git
cd ASOstool-v2
```

#### 2. Build the Docker image 

```
docker build -t asostoolv2 .
```

#### 3. Create a persistent volume for RStudio home

```
docker volume create rstudio-home
```

#### 4. Start the development container

Mount:

-   Project folder → RStudio
-   Project folder → Shiny
-   Persistent RStudio home folder

On Linux/Mac:
```         
docker run -d \
  -p 8787:8787 \
  -p 3838:3838 \
  -v $(pwd):/home/rstudio/ASOstool-v2 \
  -v $(pwd):/srv/shiny-server/ASOstool-v2 \
  -v rstudio-home:/home/rstudio \
  asostoolv2
```

On Windows (Powershell):
```
docker run -d `
  -p 8787:8787 `
  -p 3838:3838 `
  -v ${PWD}:/home/rstudio/ASOstool-v2 `
  -v ${PWD}:/srv/shiny-server/ASOstool-v2 `
  asostoolv2
```

#### 5. Log into RStudio

Open:

👉 <http://localhost:8787>

Credentials:

-   Username: rstudio
-   Password: rstudio

#### 6. Create a new R-studio project

**Create a project** in R-studio > By directory > Select ASOstool-v2 directory.

#### 7. Set up SSH keys (first time only)

Open the Terminal inside RStudio.

##### 7.1 Start SSH agent

```eval "$(ssh-agent -s)"```

##### 7.2 Generate SSH key

Replace example email with github email

```ssh-keygen -t ed25519 -C "your_email@example.com"```

Press ENTER for all prompts.

##### 7.3 Add the key to ssh-agent

```ssh-add ~/.ssh/id_ed25519```

##### 7.4 Add your public key to GitHub

```cat ~/.ssh/id_ed25519.pub```

Copy → GitHub →

**Settings** → **SSH and GPG Keys** → **New SSH key**

#### 8. Configure Git user (first time only)

```         
git config --global user.name "Your Name"`
git config --global user.email "you@example.com"
```

These persist thanks to the rstudio-home volume.

#### 9. Fix Push/Pull — convert HTTPS to SSH

If your remote url is set trough https you need to switch to SSH

```
cd /home/rstudio/ASOstool-v2
git remote set-url origin git@github.com:JimReijnhoudt/ASOstool-v2.git
```

Check:

```
git remote -v
```

You should now see:

```
origin  git@github.com:JimReijnhoudt/ASOstool-v2.git (fetch)
origin  git@github.com:JimReijnhoudt/ASOstool-v2.git (push)
```

Now the Push and Pull buttons in RStudio Git tab become active.

#### 10. Fix: “All files show as modified”

Run once:

```
git config core.autocrlf input
```

#### 11. Stop the development container

Stop the container after you're done to free up ports.
List running containers:

```
docker ps
```

Stop it:

```
docker stop <container_id>
```

Remove it (optional):

```
docker rm <container_id>
```
---
### ✔ Setup Summary Table

| Action                                                     | Required Once | Required Each Restart |
| ----------------------------------------------------------- | ------------- | --------------------- |
| Create `rstudio-home` volume                                | ✔             | ❌                     |
| Generate SSH key                                            | ✔             | ❌                     |
| Add SSH key to GitHub                                       | ✔             | ❌                     |
| Set Git username/email                                      | ✔             | ❌                     |
| Convert remote to SSH                                       | ✔             | ❌                     |
| Start ssh-agent + ssh-add (If using password for SSH key)   | ❌             | ✔                     |
| Run docker container                                        | ❌             | ✔                     |
| Stop container                                              | ❌             | ✔                     |

## 📦 Dependencies
| Component                                   | Version                | Source / Notes                                                                 |
|---------------------------------------------|------------------------|-------------------------------------------------------------------------------|
| **Base Image**                              |                        |                                                                               |
| rocker/tidyverse                            | 4.5.1                  | Docker Hub                                                                    |
| **Shiny Server**                            |                        |                                                                               |
| Shiny Server                                | 1.5.23.1030            | Installed via rocker install script                                            |
| **Manually Built Software**                 |                        |                                                                               |
| ViennaRNA                                   | 2.7.0                  | Downloaded and compiled from source                                            |
| **R Packages – CRAN**                       |                        |                                                                               |
| shiny                                       | latest                 | CRAN                                                                           |
| shinythemes                                 | latest                 | CRAN                                                                           |
| tidyverse                                   | latest                 | CRAN                                                                           |
| rlang                                       | latest                 | CRAN                                                                           |
| dplyr                                       | latest                 | CRAN                                                                           |
| DT                                          | latest                 | CRAN                                                                           |
| cluster                                     | latest                 | CRAN                                                                           |
| shinyBS                                     | latest                 | CRAN                                                                           |
| shinydashboard                              | latest                 | CRAN                                                                           |
| **R Packages – Bioconductor**               |                        |                                                                               |
| BiocManager                                | latest                 | Installed from CRAN                                                            |
| GenomicFeatures                            | latest                 | Bioconductor                                                                   |
| AnnotationDbi                              | latest                 | Bioconductor                                                                   |
| BSgenome.Hsapiens.NCBI.GRCh38              | latest                 | Bioconductor                                                                   |
| biomaRt                                    | latest                 | Bioconductor                                                                   |
| Biostrings                                 | latest                 | Bioconductor                                                                   |
| txdbmaker                                  | latest                 | Bioconductor                                                                   |
| **Generated Data**                          |                        |                                                                               |
| Homo_sapiens.GRCh38.112.gtf.gz              | Release 112            | Downloaded from Ensembl                                                       |
| txdb_hsa_biomart.db                         | generated              | Stored in `/opt/ASOstool-v2/`                                                 |
