# Daylily Design

  > Daylily was designed to be a framework to easily run, and compare multi-omic tools. Rather than one pipeline, it treats each tool as a node in a DAG, abstracting each class of tool (aligner, bam sorter, duplicate marking, SNV/SV calling, and numerous QC tools) so they may run in a combinatorial fashion. This abstraction allows for clear, and well organized, benchmarking of tools in isolation and in combination with appropriate up/downstream tooling.

## Design Principles

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
  > The process which produces the input data for daylily is critically important to the output analysis quality. Further, the analysis daylily performs is not only calling variants, but also producing batch QC reports which are invaluable to monitoring lab process health over time and serving as QC and debugging information on a per-batch viewpoint.  The consumers of these analyses are not only bioinformaticans, but Scientists processing samples, automation engineers, IT personel, technicians, operations staff, among others.  Daylily organizes results in a clear and accessible fashion so as to be of utility to all consumers of these workflows.
  
### Speed and Cost Conscious (aka- be aware of your hardware)
  > Informatics tools are notoriously challenging to run efficently. Daylily began from the bottom and worked upwards, fretting about speed & cost at every step. To get the best performance out of informatics tools, it is important to understand in which compute environments these tools will be running. Questions to ask:

  - Can a tool be re-compiled (say, with the intel compiler) to better run on modern intel hardware ( supporting the AVX-512 instruction set)? Suprisingly often, the answer is yes!
    - Is disk I/O an issue?  
      - Can the underlying disk be changed? 
      - Can I/O burden be reduced (by leveraging pip'ing between processes, or process substitution)?
    - Is sufficent monitoring in place to draw conclusions about each tools performance, and how well it is matched to the compute it is running on?
    - Are job failures handled gracefully, and re-starts of crashed workflows straight forward?
      - Are jobs running quickly enough so as to de-risk using spot instances? 
    - Are the environments to run tools reliably reproducible to instantiate (using containers or well defined virtual envs)?
