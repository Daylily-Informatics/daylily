# [DAY](https://en.wikipedia.org/wiki/Margaret_Oakley_Dayhoff)![](https://placehold.co/60x35/ff03f3/fcf2fb?text=LILLY)
_named in honor of Margaret Oakley Dahoff_


<p valign="middle"><a href=http://www.workwithcolor.com/color-converter-01.htm?cp=ff8c00><img src="docs/images/0000002.png" valign="bottom" ></a></p>

# daylily is:

## A Free Multi-omics Analysis Framework 

### (in <1hr : a $5 EC2 30x genome / a slightly more sensitive 30x genome for $10) && exhaustive QC reporting, SV callers, perf monitoring, cost reporting... 

  > Daylily connects and automates all of the systems necessary to run omic analysis reproducibly, reliablay and at massive scale in the cloud. All you need is a laptop and access to an AWS console. After a [~90m installation](docc/install/video_guide.md), you will be ready to begin processing up to thousands of genomes an hour. The code is open source and free to use! I am also absoutely available for consulting services to integrate daylily, migrate pipelines into this framework, or optimize existing pipelines.

## A Managed Service

  > Daylily offers a managed service where, depending on the analysis and TAT desired, you pay a per-sample fee for results produced by data you provide.  The interface to send sample data to daylily and for daylily to return results can be as straight forward as a properly permissioned S3 bucket. Observability around sample/batch analysis progress will allow you to know when to expect results. Contact daylily@daylilyinformatics.com for further information.

# Very General Design Overview

  [daylily_birds_eye.pdf](https://github.com/Daylily-Informatics/daylily/files/11390618/daylily_birds_eye.pdf)

## 



### Use Cases: Performance, Fscores, Costs
  >  Users may choose among the pre-validated pipeline options, or may extend the framework to run custom pipelines. The Scientific workflow manager behind the scenes happens to be [snakemake](), but any workflow manager that integrates with slurm or aws batch can extend the framework.

  - Three validated pipelines are available. Two are comprised of fully open source tools, the third leverages hardware agnostic accelerated tools from Sentieon. The pipelines and average performance across the google brain 30x Novaseq fastqs for the 7 [giab]() samples are as follows:
 
 
 | Pipeline |   SNP fscore  |  INS fscore |  DEL fscore | Indel fscore |  e2e walltime |  e2e 128vcpu instance min | Avg EC2 Cost |
 | :-------------: | :-------------: | :--------------: | :-------------: | :-------------: | :--------------: | :-------------: | :-------------: |
 |   Sentieon BWA + DNAscope (BD) |  SNP fscore | INS fscore | DEL fscore | Indel fscore | e2e walltime | e2e 128vcpu instance min | Avg EC2 Cost |
 |   BWA-MEM2 + Octopus (B2O) |  SNP fscore | INS fscore | DEL fscore | Indel fscore | e2e walltime | e2e 128vcpu instance min | Avg EC2 Cost |
 |   BWA-MEM2 + Deepvariant (B2D) |  SNP fscore | INS fscore | DEL fscore | Indel fscore | e2e walltime | e2e 128vcpu instance min | Avg EC2 Cost |



<p valign="middle"><a href=http://www.workwithcolor.com/color-converter-01.htm?cp=ff8c00><img src="docs/images/000000.png" valign="bottom" ></a></p\
>

<p valign="middle"><img src="docs/images/000000.png" valign="bottom" ></p>



### Daylily Design Considerations

#### Reproducibility (Sustainability, Portability, etc)
  > This framework fully embraces [these ideals](), and by wrapping snakemake in a similarly opinionated codebase, the hope is that these omic analysis tools will be approachable by a wider audience.
    - Daylily and it's components are thoroughly [versioned](docs/more/versioning.md), allowing for complete replication of analysis environments and results.
    - IMAGE of Entire Pipeline

#### Validation Ready
  > Informed by the requirements of running a CLIA & CAP accredited diagnostic lab, daylily will automatically produce concordance reports for arbitrary samples within the bed footprint you specify. Further, [exhaustive QC metrics]() are normalized so as to be painless to work with.

#### Meaningful Quality Metrics
  - CCCCC


<p valign="middle"><a href=http://www.workwithcolor.com/color-converter-01.htm?cp=ff8c00><img src="docs/images/000000.png" valign="bottom" ></a></p\
>

<p valign="middle"><img src="docs/images/000000.png" valign="bottom" ></p>



### Daylily Framework Details
#### SNV Calling Pipelines (validated)
##### Sentieon BWA + DNAscope (BD)
  > Picture and list of tools

##### BWA-MEM2 + Octopus (B2O)
  > Picture and	list of	tools

##### BWA-MEM2 + Deepvariant (B2D)
  >Picture and	list of	tools

#### QC Tools & Reports
  >Picture and  list of tools


#### Performance Monitoring Reports
  >Picture and  list of tools

#### Observability
  >Picture and  list of tools

#### Cost Tracking and Budget Enforcement


<p valign="middle"><a href=http://www.workwithcolor.com/color-converter-01.htm?cp=ff8c00><img src="docs/images/000000.png" valign="bottom" ></a></p>

<p valign="middle"><img src="docs/images/000000.png" valign="bottom" ></p>



# Getting Started

## Consulting
  - Depending on the level of integration, daylily should be able to be deployed and running in your environment in a matter of days.  The basic integration looks approximately like the figure below. Additionally, I am available for custom pipeline development, migration and optimization. Please contact daylily@daylilyinformatics.com if consulting services are of interest.
    - IMAGE

## Straightforward [Licensing]() 
  *With* registration (email to daylily@daylilyinformatics.com):
    - free for personal and academic use.
    - nominal fee for pre-commercial use.
    - fee for commercial use.
    - appropriate citation if publications generated results using daylily.
    - The code is otherwise available under the GNPXXX liscence ([see liscence](LICENSE)).


## Installation

### Prerequisites
  - a
  - b


## Extending
  - more to come

<p valign="middle"><a href=http://www.workwithcolor.com/color-converter-01.htm?cp=ff8c00><img src="docs/images/000000.png" valign="bottom" ></a></p>

<p valign="middle"><img src="docs/images/000000.png" valign="bottom" ></p>



<p valign="middle"><a href=http://www.workwithcolor.com/color-converter-01.htm?cp=ff8c00><img src="docs/images/000000.png" valign="bottom" ></a></p>

<p valign="middle"><a href=http://www.workwithcolor.com/color-converter-01.htm?cp=ff8c00><img src="docs/images/0000002.png" valign="bottom" ></a></p>
