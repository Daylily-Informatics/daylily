# [DAY](https://en.wikipedia.org/wiki/Margaret_Oakley_Dayhoff)![](https://placehold.co/60x35/ff03f3/fcf2fb?text=LILLY)
_named in honor of Margaret Oakley Dahoff_


<p valign="middle"><a href=http://www.workwithcolor.com/color-converter-01.htm?cp=ff8c00><img src="docs/images/0000002.png" valign="bottom" ></a></p>

# daylily is:

## A Free Multi-omics Analysis Framework 

### FQ to SNVvcf for as little as $3.50 and in under 1 hour, parallelizable to thousannds of genomes an hour
#### _PLUS_ SNV/SV calling options at other sensitivities / extensive sample + batch QC reporting / performance & cost reporting + budgeting  

  > Daylily provides a single point of contact to the myriad systems which need to be orchestrated in order to run omic analysis reproducibly, reliablay and at massive scale in the cloud. **All you need is a laptop and access to an AWS console**. After a [~90m installation](docc/install/video_guide.md), you will be ready to begin processing up to thousands of genomes an hour. 
  > This code is open source and free to use(excepting the fastest pipline)! I hope some neat tricks I deploy are of help to others [see blog](https://daylily-informatics.github.io/). 
  > [Daylily Informatics](http://daylilyinformatics.com/) is available for consulting services to integrate daylily into your operations, migrate pipelines into this framework, or optimize existing pipelines.

## A Managed Service

  > Daylily Informatics offers a managed service where, depending on the analyses and TAT desired, you pay a per-sample fee for results produced from the data you submit.  The interface to send and receive data can be as straight forward as a properly permissioned S3 bucket where new fastqs are picked up, and soon afterwards, results are returned to the same bucket. Observability around sample/batch analysis progress will allow you to know the status of your samples. Please contact [daylily@daylilyinformatics.com](https://us21.list-manage.com/contact-form?u=434d42174af0051b1571c6dce&form_id=23d28c274008c0829e07aff8d5ea2e91) for further information.

<p valign="middle"><a href=http://www.workwithcolor.com/color-converter-01.htm?cp=ff8c00><img src="docs/images/0000002.png" valign="bottom" ></a></p>

# Very General Components Overview

  > Before even getting into the cool informatics business going on, there is a boatload of complex ops systems running to manage EC2 instances, navigate spot markets, as well as mechanisms to monitor and observe all aspects of this framework. [AWS ParallelCluster](https://docs.aws.amazon.com/parallelcluster/latest/ug/what-is-aws-parallelcluster.html) is the glue holding everything together, and deserve special thanks.

  [daylily_birds_eye.pdf](https://github.com/Daylily-Informatics/daylily/files/11390618/daylily_birds_eye.pdf)
  

## Some Bioinformatics Bits, Big Picture

### The DAG For 1 Sample Running Through The `BWA-MEM2ert+Doppelmark+Deepvariant+Manta+TIDDIT+Dysgu+Svaba+QCforDays` Pipeline

   ![](docs/images/assets/ks_rg.png)
   
   - The above is actually a compressed view of the jobs managed for a sample moving through this pipeline. This view is of the dag which properly reflects parallelized jobs.... think UFO.
   
     ![](docs/images/assets/ks_dag.png)
   

#### Final Batch QC Report Produced 
  > this report is generated running the google-brain Novaseq 30x HG002 fastqs, and again downsampling to: 25,20,15,10,5x.
   [MQC](docs/images/assets/MQC_example.html)

### Performance Tuning
  > A significant amount of time was spent exploring tools, finding which environments they operated best within, and so on. I plan to use the [daylily blog]([)](https://daylily-informatics.github.io/) to discuss some of the neat insigts I was granted during this process.


## [Daylily Design Principles](docs/ops/design_principles.md)
  > Daylily was built while drawing on over 20 years of experience in clinical genomics and informatics. [These](docs/ops/design_principles.md) principles front and center while building this framework.


## Some Bioinformatics Bits, Brass Tacks


### 3 Pipelines: Performance, Fscores, Costs
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
  >x

<p valign="middle"><a href=http://www.workwithcolor.com/color-converter-01.htm?cp=ff8c00><img src="docs/images/000000.png" valign="bottom" ></a></p>

<p valign="middle"><img src="docs/images/000000.png" valign="bottom" ></p>




<p valign="middle"><a href=http://www.workwithcolor.com/color-converter-01.htm?cp=ff8c00><img src="docs/images/000000.png" valign="bottom" ></a></p>

<p valign="middle"><img src="docs/images/000000.png" valign="bottom" ></p>



<p valign="middle"><a href=http://www.workwithcolor.com/color-converter-01.htm?cp=ff8c00><img src="docs/images/000000.png" valign="bottom" ></a></p>

<p valign="middle"><a href=http://www.workwithcolor.com/color-converter-01.htm?cp=ff8c00><img src="docs/images/0000002.png" valign="bottom" ></a></p>
