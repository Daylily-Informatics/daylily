FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="b7bb3cc1ec9b03c806c2afc44963d48e3eca776fbb4eace4d38cef5e376cc774"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/alignstats_v0.2.yaml
#   prefix: /conda-envs/6f2307ffa1e0ec751250f3c8957bb9d0
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - alignstats==0.10
#     - parallel=20210822
#     - samtools
#     - pigz
#     - bcftools
#     - isa-l
#     - htslib
#     - bedtools
#     - perl
#     - tabix
RUN mkdir -p /conda-envs/6f2307ffa1e0ec751250f3c8957bb9d0
COPY workflow/envs/alignstats_v0.2.yaml /conda-envs/6f2307ffa1e0ec751250f3c8957bb9d0/environment.yaml

# Conda environment:
#   source: workflow/envs/bwasamtools_v0.1.yaml
#   prefix: /conda-envs/21acb778ff57c6be3841202554a9b391
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - r
#     - defaults
#   dependencies:
#     - samtools>=1.11
#     - bwa-mem2
#     - isa-l
#     - seqkit
#     - mbuffer
#     - jemalloc
RUN mkdir -p /conda-envs/21acb778ff57c6be3841202554a9b391
COPY workflow/envs/bwasamtools_v0.1.yaml /conda-envs/21acb778ff57c6be3841202554a9b391/environment.yaml

# Conda environment:
#   source: workflow/envs/fastp_v0.1.yaml
#   prefix: /conda-envs/8cdae75f957fca692b0938ff86cecd3f
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - r
#     - defaults
#   dependencies:
#     - fastp
#     - glances
#     - glances
#     - pigz
RUN mkdir -p /conda-envs/8cdae75f957fca692b0938ff86cecd3f
COPY workflow/envs/fastp_v0.1.yaml /conda-envs/8cdae75f957fca692b0938ff86cecd3f/environment.yaml

# Conda environment:
#   source: workflow/envs/fastqc_v0.1.yaml
#   prefix: /conda-envs/2fc8add091e996fd87194baae8c06315
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - r
#     - defaults
#   dependencies:
#     - fastqc=0.11.9
#     - pigz
#     - seqfu
#     - seqkit
RUN mkdir -p /conda-envs/2fc8add091e996fd87194baae8c06315
COPY workflow/envs/fastqc_v0.1.yaml /conda-envs/2fc8add091e996fd87194baae8c06315/environment.yaml

# Conda environment:
#   source: workflow/envs/go_left_v0.1.yaml
#   prefix: /conda-envs/46dadb26e48dcc73dc177383a287149a
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - r
#     - defaults
#   dependencies:
#     - goleft=0.2.4
RUN mkdir -p /conda-envs/46dadb26e48dcc73dc177383a287149a
COPY workflow/envs/go_left_v0.1.yaml /conda-envs/46dadb26e48dcc73dc177383a287149a/environment.yaml

# Conda environment:
#   source: workflow/envs/manta_v0.1.yaml
#   prefix: /conda-envs/41e5904c2f95d7c842003639719c6a0d
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - r
#     - defaults
#   dependencies:
#     - manta
RUN mkdir -p /conda-envs/41e5904c2f95d7c842003639719c6a0d
COPY workflow/envs/manta_v0.1.yaml /conda-envs/41e5904c2f95d7c842003639719c6a0d/environment.yaml

# Conda environment:
#   source: workflow/envs/mosdepth_v0.1.yaml
#   prefix: /conda-envs/dd55706e9804723885bc015e33761852
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - r
#     - defaults
#   dependencies:
#     - mosdepth=0.3.2
#     - r-base
#     - r-r.utils
#     - r-data.table
RUN mkdir -p /conda-envs/dd55706e9804723885bc015e33761852
COPY workflow/envs/mosdepth_v0.1.yaml /conda-envs/dd55706e9804723885bc015e33761852/environment.yaml

# Conda environment:
#   source: workflow/envs/multiqc_v0.1.yaml
#   prefix: /conda-envs/eb42f412a9d3e69ca6ddd447d25c5268
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - r
#     - defaults
#   dependencies:
#     - multiqc
RUN mkdir -p /conda-envs/eb42f412a9d3e69ca6ddd447d25c5268
COPY workflow/envs/multiqc_v0.1.yaml /conda-envs/eb42f412a9d3e69ca6ddd447d25c5268/environment.yaml

# Conda environment:
#   source: workflow/envs/peddy_v0.1.yaml
#   prefix: /conda-envs/0350a8e5d8fe59e24ac5e46f3ebf823d
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - r
#     - defaults
#   dependencies:
#     - peddy
#     - python
#     - glances
RUN mkdir -p /conda-envs/0350a8e5d8fe59e24ac5e46f3ebf823d
COPY workflow/envs/peddy_v0.1.yaml /conda-envs/0350a8e5d8fe59e24ac5e46f3ebf823d/environment.yaml

# Conda environment:
#   source: workflow/envs/picard_v0.1.yaml
#   prefix: /conda-envs/7836b4a9564d0756271066f20bab1dc3
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - r
#     - defaults
#   dependencies:
#     - picard
#     - r-base
RUN mkdir -p /conda-envs/7836b4a9564d0756271066f20bab1dc3
COPY workflow/envs/picard_v0.1.yaml /conda-envs/7836b4a9564d0756271066f20bab1dc3/environment.yaml

# Conda environment:
#   source: workflow/envs/qualimap_v0.1.yaml
#   prefix: /conda-envs/555729077ebbda9576a3ffa89971155c
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - r
#     - defaults
#   dependencies:
#     - qualimap=2.2.2d
#     - glances
RUN mkdir -p /conda-envs/555729077ebbda9576a3ffa89971155c
COPY workflow/envs/qualimap_v0.1.yaml /conda-envs/555729077ebbda9576a3ffa89971155c/environment.yaml

# Conda environment:
#   source: workflow/envs/rtgtools_v0.1.yaml
#   prefix: /conda-envs/0604b80554cb0249d7d4d126b3074466
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - r
#     - defaults
#   dependencies:
#     - rtg-tools=3.12.1
#     - parallel
#     - bcftools
#     - python=3.7
#     - pip
#     - pandas
#     - htslib
#     - fd-find
#     - pip:
#         - pyvcf3
RUN mkdir -p /conda-envs/0604b80554cb0249d7d4d126b3074466
COPY workflow/envs/rtgtools_v0.1.yaml /conda-envs/0604b80554cb0249d7d4d126b3074466/environment.yaml

# Conda environment:
#   source: workflow/envs/samtools_v0.1.yaml
#   prefix: /conda-envs/77e618d7c1ff65c0ae5625b839426e24
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - r
#     - defaults
#   dependencies:
#     - samtools
RUN mkdir -p /conda-envs/77e618d7c1ff65c0ae5625b839426e24
COPY workflow/envs/samtools_v0.1.yaml /conda-envs/77e618d7c1ff65c0ae5625b839426e24/environment.yaml

# Conda environment:
#   source: workflow/envs/sentD_v0.2.yaml
#   prefix: /conda-envs/63637cb72f6c519383211ea895252a8e
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - python==3.10
#     - samtools
#     - jemalloc
#     - htslib
#     - pigz
#     - ipython
#     - pip  # Specify pip as a dependency first
#     - pip:
#         - https://github.com/sentieon/sentieon-cli/releases/download/v1.1.0/sentieon_cli-1.1.0.tar.gz
RUN mkdir -p /conda-envs/63637cb72f6c519383211ea895252a8e
COPY workflow/envs/sentD_v0.2.yaml /conda-envs/63637cb72f6c519383211ea895252a8e/environment.yaml

# Conda environment:
#   source: workflow/envs/sentieon_v0.1.yaml
#   prefix: /conda-envs/66b4597ee54a1fae7636e766fa95d360
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - r
#     - defaults
#   dependencies:
#     - sentieon==202112.06
#     - mbuffer
#     - curl
#     - samtools
#     - jemalloc
#     - htslib
#     - pigz
#     - seqkit
#     - seqfu
#     - sambamba
#     - isa-l
RUN mkdir -p /conda-envs/66b4597ee54a1fae7636e766fa95d360
COPY workflow/envs/sentieon_v0.1.yaml /conda-envs/66b4597ee54a1fae7636e766fa95d360/environment.yaml

# Conda environment:
#   source: workflow/envs/strobe_aligner.yaml
#   prefix: /conda-envs/2c02c1c05d2cb13c29f5601fbedefa38
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - samtools
#     - isa-l
#     - seqkit
#     - mbuffer
#     - jemalloc
RUN mkdir -p /conda-envs/2c02c1c05d2cb13c29f5601fbedefa38
COPY workflow/envs/strobe_aligner.yaml /conda-envs/2c02c1c05d2cb13c29f5601fbedefa38/environment.yaml

# Conda environment:
#   source: workflow/envs/vanilla_v0.1.yaml
#   prefix: /conda-envs/1ced66bd7f328feee804fa7d9c467b48
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - r
#     - defaults
#   dependencies:
#     - parallel=20210822
#     - samtools=1.21
#     - pigz
#     - bcftools
#     - isa-l
#     - htslib
#     - bedtools
#     - perl
#     - tabix
RUN mkdir -p /conda-envs/1ced66bd7f328feee804fa7d9c467b48
COPY workflow/envs/vanilla_v0.1.yaml /conda-envs/1ced66bd7f328feee804fa7d9c467b48/environment.yaml

# Conda environment:
#   source: workflow/envs/verifybamid2_v0.1.yaml
#   prefix: /conda-envs/85335de59464b4d0307620d1b2cd4487
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - r
#     - defaults
#   dependencies:
#     - verifybamid2=2.0.1
RUN mkdir -p /conda-envs/85335de59464b4d0307620d1b2cd4487
COPY workflow/envs/verifybamid2_v0.1.yaml /conda-envs/85335de59464b4d0307620d1b2cd4487/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/6f2307ffa1e0ec751250f3c8957bb9d0 --file /conda-envs/6f2307ffa1e0ec751250f3c8957bb9d0/environment.yaml && \
    mamba env create --prefix /conda-envs/21acb778ff57c6be3841202554a9b391 --file /conda-envs/21acb778ff57c6be3841202554a9b391/environment.yaml && \
    mamba env create --prefix /conda-envs/8cdae75f957fca692b0938ff86cecd3f --file /conda-envs/8cdae75f957fca692b0938ff86cecd3f/environment.yaml && \
    mamba env create --prefix /conda-envs/2fc8add091e996fd87194baae8c06315 --file /conda-envs/2fc8add091e996fd87194baae8c06315/environment.yaml && \
    mamba env create --prefix /conda-envs/46dadb26e48dcc73dc177383a287149a --file /conda-envs/46dadb26e48dcc73dc177383a287149a/environment.yaml && \
    mamba env create --prefix /conda-envs/41e5904c2f95d7c842003639719c6a0d --file /conda-envs/41e5904c2f95d7c842003639719c6a0d/environment.yaml && \
    mamba env create --prefix /conda-envs/dd55706e9804723885bc015e33761852 --file /conda-envs/dd55706e9804723885bc015e33761852/environment.yaml && \
    mamba env create --prefix /conda-envs/eb42f412a9d3e69ca6ddd447d25c5268 --file /conda-envs/eb42f412a9d3e69ca6ddd447d25c5268/environment.yaml && \
    mamba env create --prefix /conda-envs/0350a8e5d8fe59e24ac5e46f3ebf823d --file /conda-envs/0350a8e5d8fe59e24ac5e46f3ebf823d/environment.yaml && \
    mamba env create --prefix /conda-envs/7836b4a9564d0756271066f20bab1dc3 --file /conda-envs/7836b4a9564d0756271066f20bab1dc3/environment.yaml && \
    mamba env create --prefix /conda-envs/555729077ebbda9576a3ffa89971155c --file /conda-envs/555729077ebbda9576a3ffa89971155c/environment.yaml && \
    mamba env create --prefix /conda-envs/0604b80554cb0249d7d4d126b3074466 --file /conda-envs/0604b80554cb0249d7d4d126b3074466/environment.yaml && \
    mamba env create --prefix /conda-envs/77e618d7c1ff65c0ae5625b839426e24 --file /conda-envs/77e618d7c1ff65c0ae5625b839426e24/environment.yaml && \
    mamba env create --prefix /conda-envs/63637cb72f6c519383211ea895252a8e --file /conda-envs/63637cb72f6c519383211ea895252a8e/environment.yaml && \
    mamba env create --prefix /conda-envs/66b4597ee54a1fae7636e766fa95d360 --file /conda-envs/66b4597ee54a1fae7636e766fa95d360/environment.yaml && \
    mamba env create --prefix /conda-envs/2c02c1c05d2cb13c29f5601fbedefa38 --file /conda-envs/2c02c1c05d2cb13c29f5601fbedefa38/environment.yaml && \
    mamba env create --prefix /conda-envs/1ced66bd7f328feee804fa7d9c467b48 --file /conda-envs/1ced66bd7f328feee804fa7d9c467b48/environment.yaml && \
    mamba env create --prefix /conda-envs/85335de59464b4d0307620d1b2cd4487 --file /conda-envs/85335de59464b4d0307620d1b2cd4487/environment.yaml && \
    mamba clean --all -y
