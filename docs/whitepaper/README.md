


---
---
---
# OLD NOTES -- [SEE NEW WHITEPAPER DRAFT HERE](https://github.com/Daylily-Informatics/daylily_giab_analyses/blob/main/docs/whitepaper_draft.md)
---
---
---
# [Daylily](https://github.com/Daylily-Informatics/daylily) Whitepaper

Intended to be published in f1000-research.

_analysis results being processed_

# Sections

## Abstract

## Core Topics

- Reproducibility requires repeatable compute performance.
- Cost, as driven by compute performance (and storage and network) drive decisions which influence accuracy.
- Cost should be reliably predictable (and daylily includes a comprehensive cost prediction UI and other tools: for not just compute, but data transfer, storage, and long term data management)
- Daylily offers a mechanism to share a common, and reasonably accassible to all, hardware environment & be able to better benchmark.
- WGS analysis should only cost ~ $3-4 per sample, everything daylily is doing was possible *4* years ago.
- Daylily was initially designed for tool interaction benchmarking (so, all aligners by all dedupers by all variant callers by all sv callers), and it still serves this purpose admirably.
- Compute hardware reporducibility is MUCH more of a factor when choosing to run pipelines vs. the workflow orchestrators. Yet, the orchestrators is where most attention goes.  Why? Distraction, etc.
- Workflows are largely working via snakemake. *HOWEVER* the compute framework us slurm based, and any orchestrator which can use slurm can run just as easily. *Cromwell* is working in a very rudimentary fashion, and was quick to get working. Nextflow, etc can be plugged in with minimal fuss.
- When presenting information about analysis pipelines, there should be full transparency on:
-   - user run time
    - accuracy of final results
    - expectations on vcpu time per task
    -  expectations on mem required per task
    -  COST of compute, data transfer
    -  data storage
    -  cost and availability of reference and other resources needed to run compute
    -  long term sustainability/reproducibility of tools involved
    -   full details on tool versions, etc,...
 
## Daylily Design
This is going to end up 2 repos, one for framework stuff, another for analysis workflows.

### Infrasttructure Framework

### Analysis Suite

### Integrated Cost / Budget Prediction, Monitoring, Safeguards, Process Interactions

## Daylily Performance


### Dataset
- All 7 GIAB samples, google brain 30x 2x150 ILMN no amp WGS fastqs.
- `b37` human reference
  - *105* SV vcfs produced for analysis (7GIAB * 3 aligners * 1 deduper * 5 SV callers).
  - *63* SNV vcfs produced for analysis  (7GIAB * 3 aligners * 1 deduper * 3 SV callers)..
  - *A LOT* of qc reports, rolled into *1* multiqc report.

- `hg38` human reference
  - *105* SV vcfs produced for analysis (7GIAB * 3 aligners * 1 deduper * 5 SV callers).
  - *63* SNV vcfs produced for analysis  (7GIAB * 3 aligners * 1 deduper * 3 SV callers)..
  - *A LOT* of qc reports, rolled into *1* multiqc report.


### Cost & Accuracy Performance
- *~ $3-4* per 30x genome, in ~90min walltime. NOTE: there are still opportunities to reduce this by ~20%, given time/resources. NOTE2: raw compute cost can be reduced by ~$1 using sentieon, however, there will be additional s/w lisc fees to add to this.
- Fscores acheivable are as good or better than expectations from publications in the literature (so `0.999+` ), depending on tool combination.

### Toolsets
All of these have been optimized for the daylily compute framework instances.

#### Aligners
- bwa mem2
- strobealigner
- sentieon bwa mem

#### Deduplication
- doppelmark

#### SNV Callers
- deepvariant
- clair3
- lofreq2
- octopus
- sentieon DNAscope

##### Automated SNV Concordance Reports
- for all SNV vcfs, there is a detailed analysis of concordance.

#### SV Callers
- manta
- tiddit
- dysgu
- svaba _being tested still_

##### Automated SV Concordance Repors
_to be done... assistance welcome~_ 


### QC Reporting
A large number of QC tools are enabled by default, and aggregated into a final *multiqc* report
- peddy
- qualimap
- fastqc
- cov eveness
- mosdepth
- goleft
- verifybam
- alignstats
- vcfstats
- bamstats




