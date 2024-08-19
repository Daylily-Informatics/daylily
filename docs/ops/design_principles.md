# Daylily Design

  > Daylily was designed to be a framework to easily run, and compare, multi-omic tools. Rather than one pipeline, it treats each tool as a node in a DAG, abstracting each class of tool (aligner, bam sorter, duplicate marking, SNV/SV calling, and numerous QC tools) so they may run in a combinatorial fashion. This abstraction allows for clear and well organized benchmarking of tools in isolation and in combination with compatible up/downstream tooling. Once a set of tools is found useful as a _pipeline_, it can be formalized and all compute resources and software is fully versioned. This detailed versioning will satisfy internal quality and external regulatory requirements. Automated concordance reporting further streamlines validations of your chosen pipeline.

## Principles

### [Sustainable](https://f1000research.com/articles/10-33/v2)
Defined as:
  - Reproducible
  - Automated
  - Scalable
  - Portable
  - Transparent
  - Adaptable
  - Readable
  - Traceable
  - Documented
  
### Process & End User Aware
  > The analysis daylily performs is not only calling variants, but also producing batch QC reports which are invaluable to monitoring lab process health over time and providing QC and debugging information on a per-batch and per-sample basis.  The consumers of these analyses are not only bioinformaticians, but Scientists processing samples, automation engineers, IT personnel, technicians, operations staff, among others.  Daylily organizes results in a clear and accessible fashion so as to be of utility to all consumers of these workflows outputs.
  
### Usability, Speed & Cost Conscious
  > Informatics tools are notoriously challenging to run efficiently. Daylily began from the bottom and worked upwards, fretting about speed & cost at every step. To get the best performance out of informatics tools, it is important to understand the compute environments these tools will be running. Questions to ask:

  - Can a tool be re-compiled (say, with the intel compiler) to better run on modern intel hardware ( supporting the AVX-512 instruction set)? Surprisingly often, the answer is yes!
  - Is disk I/O an issue?  
    - Can the underlying disk be changed?
    - Can I/O burden be reduced (by leveraging pip'ing between processes, or process substitution)?
  - Is sufficient monitoring in place to draw conclusions about each tool's performance? Are these data produced in-line with all analyses, and accessible to any user to review?
  - Are job failures handled gracefully, and re-starts of failed workflows straight forward?
    - Are jobs running quickly enough so as to de-risk using spot instances? 
  - Are the environments to run tools reliably reproducible to instantiate (using containers or well defined virtual envs)?
  - Are technology choices made with preference given to solutions which are widely understood and not boutique?
  - Can 1,000's of samples be run in parallel? (theoretically, yes!)
  - Can this code be reasonably expected to still run in 1, 5, 10 years? ( packaging and well defined hardware specs are needed )
  - Is there sufficient logging and observability to debug problems when they arise?
  - Can costs be allocated by project to all resources required to produce these analyses? Are cost reports easy to produce? Can budgets be set by project?
  
### Clinical Applications In Mind
  - Daylily has been built with the needs of clinical users in mind. To this end, concordance analysis is automated and flexible. 
    - Users may specify, per-sample, a bed region to analyze and truth vcf & a concordance report ( broken down by: SNPts, SNPtv, INS, DEL, Indel) will automatically be produced. Further, the pipelines daily has formalized all come with full analysis results, concordance reports, and the tools to reproduce these results in your own hands.
  - Genome builds `b37`, `hg38` are both supported. `T2T` is in testing.
  
### Managed Analysis, or Run Yourself
  > Daylily can offer a managed analysis service, or license the code for use in your own compute environment ( free with registration at daylily@daylilyinformatics.com for personal/academic use, for commercial uses contact daylily@daylilyinformatics.com ).


