# [DAYLILY](http://daylilyinformatics.com/)

## Free, Fast(~60m) & Frugal($3.5 EC2) Multi-omics Analysis
<p valign="middle"><a href=http://www.workwithcolor.com/color-converter-01.htm?cp=ff8c00><img src="docs/images/0000002.png" valign="bottom" ></a></p>
30x FQ to SNVvcf at $3.50 EC2 costs, completes in ~60m & process thousannds of genomes an hour

- _PLUS_ SNV/SV calling options at other sensitivities / extensive sample + batch QC reporting / performance & cost reporting + budgeting  

- Daylily provides a single point of contact to the myriad systems which need to be orchestrated in order to run omic analysis reproducibly, reliablay and at massive scale in the cloud. **All you need is a laptop and access to an AWS console**. After a [~90m installation](docs/install/video_guide.md), you will be ready to begin processing up to thousands of genomes an hour. 

- Daylily is open source and free to use(excepting the fastest pipline)! I hope some neat tricks I deploy are of help to others [see blog](https://daylily-informatics.github.io/). 

  > **Note**
  > [Daylily Informatics](http://daylilyinformatics.com/) is available for consulting services to integrate daylily into your operations, migrate pipelines into this framework, or optimize existing pipelines. [daylily@daylilyinformatics.com](https://us21.list-manage.com/contact-form?u=434d42174af0051b1571c6dce&form_id=23d28c274008c0829e07aff8d5ea2e91)

## Managed Analysis Service
- Daylily Informatics offers a managed genomic analysis service where, depending on the analyses and TAT desired, you pay a per-sample fee for daylily to run the desired analysis.
- The gist of the standard deployment can be reviewed [here](https://github.com/Daylily-Informatics/daylily/blob/main/README.md#managed-genomics-analysis-services). 
  
- Please contact [daylily@daylilyinformatics.com](https://us21.list-manage.com/contact-form?u=434d42174af0051b1571c6dce&form_id=23d28c274008c0829e07aff8d5ea2e91) for further information.


  <p valign="middle"><a href=http://www.workwithcolor.com/color-converter-01.htm?cp=ff8c00><img src="docs/images/000000.png" valign="bottom" ></a></p>

  <p valign="middle"><img src="docs/images/000000.png" valign="bottom" ></p>

# General Components Overview

  > Before even getting into the cool informatics business going on, there is a boatload of complex ops systems running to manage EC2 instances, navigate spot markets, as well as mechanisms to monitor and observe all aspects of this framework. [AWS ParallelCluster](https://docs.aws.amazon.com/parallelcluster/latest/ug/what-is-aws-parallelcluster.html) is the glue holding everything together, and deserve special thanks.
  
![DEC_components_v2](https://user-images.githubusercontent.com/4713659/236144817-d9b26d68-f50b-423b-8e46-410b05911b12.png)

# Managed Genomics Analysis Services
The system is designed to be robust, secure, auditable, and should only take a matter of days to stand up. [Please contact me for further details](https://us21.list-manage.com/contact-form?u=434d42174af0051b1571c6dce&form_id=23d28c274008c0829e07aff8d5ea2e91).


![daylily_managed_service](https://user-images.githubusercontent.com/4713659/236186668-6ea2ec81-9fe4-4549-8ed0-6fcbd4256dd4.png)


<p valign="middle"><img src="docs/images/000000.png" valign="bottom" ></p>

## Some Bioinformatics Bits, Big Picture

### The DAG For 1 Sample Running Through The `BWA-MEM2ert+Doppelmark+Deepvariant+Manta+TIDDIT+Dysgu+Svaba+QCforDays` Pipeline

   ![](docs/images/assets/ks_rg.png)
   
   - The above is actually a compressed view of the jobs managed for a sample moving through this pipeline. This view is of the dag which properly reflects parallelized jobs.... think UFO.
   
     ![](docs/images/assets/ks_dag.png)


## [Daylily Design Principles](docs/ops/design_principles.md)
  > Daylily was built while drawing on over 20 years of experience in clinical genomics and informatics. [These](docs/ops/design_principles.md) principles were kept front and center while building this framework.


## Some Bioinformatics Bits, Brass Tacks

### Three Pipelines: Performance, Fscores, Costs
  >  Users may choose among the pre-validated pipeline options, or may extend the framework to run custom pipelines. The Scientific workflow manager behind the scenes happens to be [snakemake](), but any workflow manager that integrates with slurm or aws batch can extend the framework.

  - Three validated pipelines are available. Two are comprised of fully open source tools, the third leverages hardware agnostic accelerated tools from Sentieon. The pipelines and average performance across the google brain 30x Novaseq fastqs for the 7 [giab]() samples are as follows:
 
 
 | Pipeline |   SNP fscore  |  INS fscore |  DEL fscore | Indel fscore |  e2e walltime |  e2e 128vcpu instance min | Avg EC2 Cost |
 | :-------------: | :-------------: | :--------------: | :-------------: | :-------------: | :--------------: | :-------------: | :-------------: |
 |   Sentieon BWA + DNAscope (BD) |  SNP fscore | INS fscore | DEL fscore | Indel fscore | e2e walltime | e2e 128vcpu instance min | Avg EC2 Cost |
 |   BWA-MEM2 + Octopus (B2O) |  SNP fscore | INS fscore | DEL fscore | Indel fscore | e2e walltime | e2e 128vcpu instance min | Avg EC2 Cost |
 |   BWA-MEM2 + Deepvariant (B2D) |  SNP fscore | INS fscore | DEL fscore | Indel fscore | e2e walltime | e2e 128vcpu instance min | Avg EC2 Cost |




<p valign="middle"><a href=http://www.workwithcolor.com/color-converter-01.htm?cp=ff8c00><img src="docs/images/000000.png" valign="bottom" ></a></p\>

<p valign="middle"><img src="docs/images/000000.png" valign="bottom" ></p>



### Daylily Framework, More

#### [Batch QC HTML Summary Report](http://daylilyinformatics.com:8081/components/daylily_qc_reports/reports/DAY_final_multiqc.html)
> Batch comprised of google-brain Novaseq 30x HG002 fastqs, and again downsampling to: 25,20,15,10,5x.     

![](docs/images/assets/day_qc_1.png)

![](docs/images/assets/day_qc_2.png)
    
    
### [Consistent + Easy To Navigate Results Directory & File Structure](/docs/ops/dir_and_file_scheme.md)
   - A visualization of just the directories (minus log dirs) created by daylily
    _b37 shown, hg38 is supported as well_
    [![](docs/images/assets/datlily_tree.png)](docs/ops/tree.md)
     - [with files](docs/ops/tree_full.png)     
    
### [Automated Concordance Analysis Table](http://daylilyinformatics.com:8081/components/daylily_qc_reports/other_reports/giabhcr_concordance_mqc.tsv)
  > Reported fasceted by: SNPts, SNPtv, INS>0-<51, DEL>0-51, Indel>0-<51.
  > Generated when the correct info is set in the analysis_manifest.


#### [Performance Monitoring Reports]()
  > Picture and  list of tools

#### [Observability w/CloudWatch Dashboard)[https://us-east-2.console.aws.amazon.com/cloudwatch/home?region=us-east-2#]
  > ![](docs/images/assets/cloudwatch.png)
  > ![](/docs/images/assets.cloudwatch_2.png)
  > ![](/docs/images/assets.cloudwatch3.png)

#### [Cost Tracking and Budget Enforcement])(https://aws.amazon.com/blogs/compute/using-cost-allocation-tags-with-aws-parallelcluster/)
  > ![](https://d2908q01vomqb2.cloudfront.net/1b6453892473a467d07372d45eb05abc2031647a/2020/07/23/Billing-console-projects-grouping.png)
  - ![](docs/images/assets/costs1.png)
  - ![](docs/images/assets/costs2.png)
  
  
<p valign="middle"><a href=http://www.workwithcolor.com/color-converter-01.htm?cp=ff8c00><img src="docs/images/000000.png" valign="bottom" ></a></p>

# Coming Sooner Than Latter
- Structural Variant Calling Concordance Analysis For The SV Callers:
  - Manta
  - TIDDIT
  - Svaba
  - Dysgu
  - Octopus (which is a fantastic small SV caller!)
- Annotation of SNV / SV `vcf` files with potentially clinically relevant info (VEP is in testing).
- Document the steps to quickly re-run the 7 30x GIAB samples from scratch.
- Explore hybrid assemblies using short and long reads (ONT + PacBio).


# DOCUMENTATION (WIP)

## [Installation](docs/install/Install.md)

## [Cost Tagging](docs/ops/cost_tagging.md)

## [DEC Config](docs/ops/config.md)

## [dy-CLI](docs/ops/dycli.md)

## [Analysis Manifest](docs/ops/analysis_manifest.md)

## [Batch Quality Control](docs/ops/batch_quality_control.md)

## [Visualizations](docs/ops/visualizations.md)

## [Running Tests](docs/ops/tests.md)

<p valign="middle"><img src="docs/images/000000.png" valign="bottom" ></p>

<p valign="middle"><a href=http://www.workwithcolor.com/color-converter-01.htm?cp=ff8c00><img src="docs/images/000000.png" valign="bottom" ></a></p>

<p valign="middle"><img src="docs/images/000000.png" valign="bottom" ></p>



<p valign="middle"><a href=http://www.workwithcolor.com/color-converter-01.htm?cp=ff8c00><img src="docs/images/000000.png" valign="bottom" ></a></p>

<p valign="middle"><a href=http://www.workwithcolor.com/color-converter-01.htm?cp=ff8c00><img src="docs/images/0000002.png" valign="bottom" ></a></p>

# [DAY](https://en.wikipedia.org/wiki/Margaret_Oakley_Dayhoff)![](https://placehold.co/60x35/ff03f3/fcf2fb?text=LILLY)

_named in honor of Margaret Oakley Dahoff_
