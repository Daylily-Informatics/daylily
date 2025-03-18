FROM continuumio/miniconda3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="82749d950ab686518a025ca0575ffbdc93398dc214d860067c8d23c9b5007b33"

RUN conda install -n base conda=25.1.1 -c conda-forge -y && \
    conda clean -afy

RUN apt-get update && apt-get install -y curl && rm -rf /var/lib/apt/lists/*

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/alignstats_v0.2.yaml
#   prefix: /conda-envs/49487071b4ddc399777cf2350c57bc4e
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
#     - isa-l==2.31.0
#     - htslib
#     - bedtools
#     - perl
#     - tabix
RUN mkdir -p /conda-envs/49487071b4ddc399777cf2350c57bc4e
COPY workflow/envs/alignstats_v0.2.yaml /conda-envs/49487071b4ddc399777cf2350c57bc4e/environment.yaml

# Conda environment:
#   source: workflow/envs/bwasamtools_v0.1.yaml
#   prefix: /conda-envs/4579c0ed9c9c4d96084f2864869f1063
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - r
#     - defaults
#   dependencies:
#     - samtools>=1.11
#     - bwa-mem2
#     - isa-l==2.31.0
#     - seqkit
#     - mbuffer
RUN mkdir -p /conda-envs/4579c0ed9c9c4d96084f2864869f1063
COPY workflow/envs/bwasamtools_v0.1.yaml /conda-envs/4579c0ed9c9c4d96084f2864869f1063/environment.yaml

# Conda environment:
#   source: workflow/envs/doppelmark_v0.1.yaml
#   prefix: /conda-envs/dd175262a8290a96c4741e28b33a9e7d
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - parallel=20210822
#     - samtools=1.21
#     - pigz
#     - bcftools
#     - isa-l==2.31.0
#     - htslib
#     - bedtools
#     - perl
#     - tabix
#     - mbuffer
RUN mkdir -p /conda-envs/dd175262a8290a96c4741e28b33a9e7d
COPY workflow/envs/doppelmark_v0.1.yaml /conda-envs/dd175262a8290a96c4741e28b33a9e7d/environment.yaml

# Conda environment:
#   source: workflow/envs/dysgu_sv_v0.2.yaml
#   prefix: /conda-envs/6e967b3f65887c18414bac9bbc2db597
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - dysgu==1.8.1
RUN mkdir -p /conda-envs/6e967b3f65887c18414bac9bbc2db597
COPY workflow/envs/dysgu_sv_v0.2.yaml /conda-envs/6e967b3f65887c18414bac9bbc2db597/environment.yaml

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
#   source: workflow/envs/lofreq2_v0.1.yaml
#   prefix: /conda-envs/ccb6b0743ce0883290a3cb2d8f2b378f
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - r
#     - defaults
#   dependencies:
#     - lofreq
#     - pigz
#     - bcftools
RUN mkdir -p /conda-envs/ccb6b0743ce0883290a3cb2d8f2b378f
COPY workflow/envs/lofreq2_v0.1.yaml /conda-envs/ccb6b0743ce0883290a3cb2d8f2b378f/environment.yaml

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
#   prefix: /conda-envs/a79b2570d62539822891587754f18437
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - multiqc==1.27.1
RUN mkdir -p /conda-envs/a79b2570d62539822891587754f18437
COPY workflow/envs/multiqc_v0.1.yaml /conda-envs/a79b2570d62539822891587754f18437/environment.yaml

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
#   prefix: /conda-envs/11b98826e45e7ea36ff5714da7e0dcc2
#   channels:
#     - conda-forge
#     - bioconda
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
RUN mkdir -p /conda-envs/11b98826e45e7ea36ff5714da7e0dcc2
COPY workflow/envs/rtgtools_v0.1.yaml /conda-envs/11b98826e45e7ea36ff5714da7e0dcc2/environment.yaml

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
#   prefix: /conda-envs/70cd4bd7c7283fa0aa63c5db17ad2443
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
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
RUN mkdir -p /conda-envs/70cd4bd7c7283fa0aa63c5db17ad2443
COPY workflow/envs/sentD_v0.2.yaml /conda-envs/70cd4bd7c7283fa0aa63c5db17ad2443/environment.yaml

# Conda environment:
#   source: workflow/envs/sentieon_v0.1.yaml
#   prefix: /conda-envs/fe530aa03f9b5e73e03498b40b1949ef
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
#     - isa-l==2.31.0
RUN mkdir -p /conda-envs/fe530aa03f9b5e73e03498b40b1949ef
COPY workflow/envs/sentieon_v0.1.yaml /conda-envs/fe530aa03f9b5e73e03498b40b1949ef/environment.yaml

# Conda environment:
#   source: workflow/envs/strobe_aligner.yaml
#   prefix: /conda-envs/8f717e76e5baba8bcba98d20b89cdb54
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - samtools
#     - isa-l==2.31.0
#     - seqkit
#     - mbuffer
RUN mkdir -p /conda-envs/8f717e76e5baba8bcba98d20b89cdb54
COPY workflow/envs/strobe_aligner.yaml /conda-envs/8f717e76e5baba8bcba98d20b89cdb54/environment.yaml

# Conda environment:
#   source: workflow/envs/vanilla_v0.1.yaml
#   prefix: /conda-envs/ec396639729521c9b983e237c639986c
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
#     - isa-l==2.31.0
#     - htslib
#     - bedtools
#     - perl
#     - tabix
RUN mkdir -p /conda-envs/ec396639729521c9b983e237c639986c
COPY workflow/envs/vanilla_v0.1.yaml /conda-envs/ec396639729521c9b983e237c639986c/environment.yaml

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

RUN conda env create --prefix /conda-envs/49487071b4ddc399777cf2350c57bc4e --file /conda-envs/49487071b4ddc399777cf2350c57bc4e/environment.yaml && \
    conda env create --prefix /conda-envs/4579c0ed9c9c4d96084f2864869f1063 --file /conda-envs/4579c0ed9c9c4d96084f2864869f1063/environment.yaml && \
    conda env create --prefix /conda-envs/dd175262a8290a96c4741e28b33a9e7d --file /conda-envs/dd175262a8290a96c4741e28b33a9e7d/environment.yaml && \
    conda env create --prefix /conda-envs/6e967b3f65887c18414bac9bbc2db597 --file /conda-envs/6e967b3f65887c18414bac9bbc2db597/environment.yaml && \
    conda env create --prefix /conda-envs/2fc8add091e996fd87194baae8c06315 --file /conda-envs/2fc8add091e996fd87194baae8c06315/environment.yaml && \
    conda env create --prefix /conda-envs/46dadb26e48dcc73dc177383a287149a --file /conda-envs/46dadb26e48dcc73dc177383a287149a/environment.yaml && \
    conda env create --prefix /conda-envs/ccb6b0743ce0883290a3cb2d8f2b378f --file /conda-envs/ccb6b0743ce0883290a3cb2d8f2b378f/environment.yaml && \
    conda env create --prefix /conda-envs/41e5904c2f95d7c842003639719c6a0d --file /conda-envs/41e5904c2f95d7c842003639719c6a0d/environment.yaml && \
    conda env create --prefix /conda-envs/dd55706e9804723885bc015e33761852 --file /conda-envs/dd55706e9804723885bc015e33761852/environment.yaml && \
    conda env create --prefix /conda-envs/a79b2570d62539822891587754f18437 --file /conda-envs/a79b2570d62539822891587754f18437/environment.yaml && \
    conda env create --prefix /conda-envs/0350a8e5d8fe59e24ac5e46f3ebf823d --file /conda-envs/0350a8e5d8fe59e24ac5e46f3ebf823d/environment.yaml && \
    conda env create --prefix /conda-envs/7836b4a9564d0756271066f20bab1dc3 --file /conda-envs/7836b4a9564d0756271066f20bab1dc3/environment.yaml && \
    conda env create --prefix /conda-envs/555729077ebbda9576a3ffa89971155c --file /conda-envs/555729077ebbda9576a3ffa89971155c/environment.yaml && \
    conda env create --prefix /conda-envs/11b98826e45e7ea36ff5714da7e0dcc2 --file /conda-envs/11b98826e45e7ea36ff5714da7e0dcc2/environment.yaml && \
    conda env create --prefix /conda-envs/77e618d7c1ff65c0ae5625b839426e24 --file /conda-envs/77e618d7c1ff65c0ae5625b839426e24/environment.yaml && \
    conda env create --prefix /conda-envs/70cd4bd7c7283fa0aa63c5db17ad2443 --file /conda-envs/70cd4bd7c7283fa0aa63c5db17ad2443/environment.yaml && \
    conda env create --prefix /conda-envs/fe530aa03f9b5e73e03498b40b1949ef --file /conda-envs/fe530aa03f9b5e73e03498b40b1949ef/environment.yaml && \
    conda env create --prefix /conda-envs/8f717e76e5baba8bcba98d20b89cdb54 --file /conda-envs/8f717e76e5baba8bcba98d20b89cdb54/environment.yaml && \
    conda env create --prefix /conda-envs/ec396639729521c9b983e237c639986c --file /conda-envs/ec396639729521c9b983e237c639986c/environment.yaml && \
    conda env create --prefix /conda-envs/85335de59464b4d0307620d1b2cd4487 --file /conda-envs/85335de59464b4d0307620d1b2cd4487/environment.yaml && \
    conda clean --all -y
