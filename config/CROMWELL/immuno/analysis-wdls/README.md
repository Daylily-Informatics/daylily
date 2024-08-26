# Workflows for Genomic Analysis

A collection of WDL workflows for analysis of genomic sequencing data. These enable primary analysis of many common experiments, including somatic and germline variant calling, RNA sequencing, epigenomic assays, and ommunogenomics approaches for cancer vaccine development.

These pipelines were developed by the [Washington University Division of Oncology](https://oncology.wustl.edu/), in collaboration with the [McDonnell Genome Institute](https://genome.wustl.edu/) and other members of the [Washington University School of Medicine](https://medicine.wustl.edu/).  Much of the code and structure is adapted from a sister repository of [workflows in the CWL language](https://github.com/genome/analysis-workflows/).


## Getting Started

### Workflows
Download our repository with `git clone https://github.com/wustl-oncology/analysis-wdls.git`.

Terra hosts an excellent guide to [getting started with WDL](https://support.terra.bio/hc/en-us/articles/360037117492-Overview-Getting-started-with-WDL).

### Workflow Execution
These are primarily designed to be run with the cromwell workflow system, either on a self-hosted server through the Google Cloud Platform ([details here](https://github.com/wustl-oncology/cloud-workflows)), through [Terra](https://terra.bio/)), or on local compute clusters. They have been extensively used in these contexts, though the portability of WDL means that running on other platforms should be straightforward as well.

### Docker
In order to provide a portable environment, each tool in our workflow has a designated Docker container. [Download Docker here](https://www.docker.com/products/docker-desktop).

### Data
Input "bundles" containing reference genomes and annotations are available for use/download along with documentation on how they were created.

### Help
We answer questions through github issues on this repository, and also have compiled a list of [common errors](https://github.com/wustl-oncology/analysis-wdls/blob/main/docs/common_errors.md) that may be useful.

## Contributions

A big thanks to all of the developers, bioinformaticians, and scientists who built this resource. For a complete list of software contributions, i.e. commits, to this repository, please see the GitHub Contributors both to [this repository](https://github.com/wustl-oncology/analysis-wdls/graphs/contributors) as well as to the [analysis-workflows](https://github.com/genome/analysis-workflows/graphs/contributors) repo.

Genomic data evolves rapidly and pull requests and improvements that add new features are welcome!
 
