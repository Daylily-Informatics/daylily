FROM condaforge/condaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="4f57563e7b653fd7ad7fcc0c65e1ce7cfb6e9b3c41c6395de03cf63291da4b01"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/bwasamtools_v0.1.yaml
#   prefix: /conda-envs/e332402e4f0181c4533d4808ed2bdab5
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - r

#     - defaults
#   dependencies:
#     - samtools>=1.11
#     - bwa-mem2
#     - pigz
#     - biobambam
#     - seqkit
#     - seqfu
#     - sambamba
RUN mkdir -p /conda-envs/e332402e4f0181c4533d4808ed2bdab5
COPY workflow/envs/bwasamtools_v0.1.yaml /conda-envs/e332402e4f0181c4533d4808ed2bdab5/environment.yaml

# Conda environment:
#   source: workflow/envs/fastp_v0.1.yaml
#   prefix: /conda-envs/2c1fcdef6fd8b78abb453cf367d33af6
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
RUN mkdir -p /conda-envs/2c1fcdef6fd8b78abb453cf367d33af6
COPY workflow/envs/fastp_v0.1.yaml /conda-envs/2c1fcdef6fd8b78abb453cf367d33af6/environment.yaml

# Conda environment:
#   source: workflow/envs/fastqc_v0.1.yaml
#   prefix: /conda-envs/1f6a7e34ff1bf45e82e319f82e1f3908
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
RUN mkdir -p /conda-envs/1f6a7e34ff1bf45e82e319f82e1f3908
COPY workflow/envs/fastqc_v0.1.yaml /conda-envs/1f6a7e34ff1bf45e82e319f82e1f3908/environment.yaml

# Conda environment:
#   source: workflow/envs/go_left_v0.1.yaml
#   prefix: /conda-envs/d9ae57b0259467e56b5ec3edd572f378
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - r

#     - defaults
#   dependencies:
#     - goleft=0.2.4
RUN mkdir -p /conda-envs/d9ae57b0259467e56b5ec3edd572f378
COPY workflow/envs/go_left_v0.1.yaml /conda-envs/d9ae57b0259467e56b5ec3edd572f378/environment.yaml

# Conda environment:
#   source: workflow/envs/hisat2_v0.1.yaml
#   prefix: /conda-envs/13d6c08bdbd8adaab139b9fd6e0ca041
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - r

#     - defaults
#   dependencies:
#       - hisat2=2.2.1
#       - bedtools
#       - htslib
#       - pigz
#       - samtools
RUN mkdir -p /conda-envs/13d6c08bdbd8adaab139b9fd6e0ca041
COPY workflow/envs/hisat2_v0.1.yaml /conda-envs/13d6c08bdbd8adaab139b9fd6e0ca041/environment.yaml

# Conda environment:
#   source: workflow/envs/mosdepth_v0.1.yaml
#   prefix: /conda-envs/eafc04c861b77d27fd56364ce234b91c
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
RUN mkdir -p /conda-envs/eafc04c861b77d27fd56364ce234b91c
COPY workflow/envs/mosdepth_v0.1.yaml /conda-envs/eafc04c861b77d27fd56364ce234b91c/environment.yaml

# Conda environment:
#   source: workflow/envs/multiqc_v0.1.yaml
#   prefix: /conda-envs/4076a60bcc18d220034a09af94d1ce9d
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - r

#     - defaults
#   dependencies:
#     - multiqc
RUN mkdir -p /conda-envs/4076a60bcc18d220034a09af94d1ce9d
COPY workflow/envs/multiqc_v0.1.yaml /conda-envs/4076a60bcc18d220034a09af94d1ce9d/environment.yaml

# Conda environment:
#   source: workflow/envs/octopus_v0.7.yaml
#   prefix: /conda-envs/a4543ea14359b1714e597b0e70c6a64b
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - r

#     - defaults
#   dependencies:
#     - _libgcc_mutex=0.1
#     - _openmp_mutex=4.5
#     - boost-cpp=1.74.0
#     - bzip2=1.0.8
#     - c-ares=1.18.1
#     - ca-certificates=2022.12.7
#     - gmp=6.2.1
#     - htslib=1.14
#     - icu=69.1
#     - keyutils=1.6.1
#     - krb5=1.20.1
#     - libcurl=7.87.0
#     - libdeflate=1.10
#     - libedit=3.1.20191231
#     - libev=4.33
#     - libgcc-ng=12.2.0
#     - libgomp=12.2.0
#     - libnghttp2=1.51.0
#     - libssh2=1.10.0
#     - libstdcxx-ng=12.2.0
#     - libzlib=1.2.13
#     - ncurses=6.3
#     - octopus=0.7.4
#     - openssl=1.1.1s
#     - xz=5.2.6
#     - zlib=1.2.13
#     - zstd=1.5.2
#     - htslib
RUN mkdir -p /conda-envs/a4543ea14359b1714e597b0e70c6a64b
COPY workflow/envs/octopus_v0.7.yaml /conda-envs/a4543ea14359b1714e597b0e70c6a64b/environment.yaml

# Conda environment:
#   source: workflow/envs/peddy_v0.1.yaml
#   prefix: /conda-envs/753ea3fe000cd80431db3f4fae0ccbe2
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
RUN mkdir -p /conda-envs/753ea3fe000cd80431db3f4fae0ccbe2
COPY workflow/envs/peddy_v0.1.yaml /conda-envs/753ea3fe000cd80431db3f4fae0ccbe2/environment.yaml

# Conda environment:
#   source: workflow/envs/picard_v0.1.yaml
#   prefix: /conda-envs/7b17fcd547d9025990813fe914db3286
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - r

#     - defaults
#   dependencies:
#     - picard
#     - r-base
RUN mkdir -p /conda-envs/7b17fcd547d9025990813fe914db3286
COPY workflow/envs/picard_v0.1.yaml /conda-envs/7b17fcd547d9025990813fe914db3286/environment.yaml

# Conda environment:
#   source: workflow/envs/qualimap_v0.1.yaml
#   prefix: /conda-envs/9bde885de7e411ded58cabb1dae00fb9
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - r

#     - defaults
#   dependencies:
#     - qualimap=2.2.2d
#     - glances
RUN mkdir -p /conda-envs/9bde885de7e411ded58cabb1dae00fb9
COPY workflow/envs/qualimap_v0.1.yaml /conda-envs/9bde885de7e411ded58cabb1dae00fb9/environment.yaml

# Conda environment:
#   source: workflow/envs/rtgtools_v0.1.yaml
#   prefix: /conda-envs/4e1fe845949f2d4e6edd63dc1c2a3793
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
RUN mkdir -p /conda-envs/4e1fe845949f2d4e6edd63dc1c2a3793
COPY workflow/envs/rtgtools_v0.1.yaml /conda-envs/4e1fe845949f2d4e6edd63dc1c2a3793/environment.yaml

# Conda environment:
#   source: workflow/envs/samtools_v0.1.yaml
#   prefix: /conda-envs/f080e4234b0c70dc690b6c52fbfe4fbe
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - r

#     - defaults
#   dependencies:
#     - samtools
RUN mkdir -p /conda-envs/f080e4234b0c70dc690b6c52fbfe4fbe
COPY workflow/envs/samtools_v0.1.yaml /conda-envs/f080e4234b0c70dc690b6c52fbfe4fbe/environment.yaml

# Conda environment:
#   source: workflow/envs/sentD_v0.1.yaml
#   prefix: /conda-envs/244b237b7cd6c4b6fef25b5fcd603b7a
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - r

#     - defaults
#   dependencies:
#     - sentieon==202112.06
#     - samtools
#     - jemalloc
#     - htslib
#     - pigz
#     - seqkit
#     - seqfu
RUN mkdir -p /conda-envs/244b237b7cd6c4b6fef25b5fcd603b7a
COPY workflow/envs/sentD_v0.1.yaml /conda-envs/244b237b7cd6c4b6fef25b5fcd603b7a/environment.yaml

# Conda environment:
#   source: workflow/envs/sentieon_v0.1.yaml
#   prefix: /conda-envs/e40c66cc926b5d83dd51ff71746f8a18
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - r

#     - defaults
#   dependencies:
#     - sentieon==202112.06
#     - samtools
#     - jemalloc
#     - htslib
#     - pigz
#     - seqkit
#     - seqfu
#     - sambamba
RUN mkdir -p /conda-envs/e40c66cc926b5d83dd51ff71746f8a18
COPY workflow/envs/sentieon_v0.1.yaml /conda-envs/e40c66cc926b5d83dd51ff71746f8a18/environment.yaml

# Conda environment:
#   source: workflow/envs/vanilla_v0.1.yaml
#   prefix: /conda-envs/4b68e41282e18d46ed0c01f28070afa8
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - r

#     - defaults
#   dependencies:
#     - parallel=20210822
#     - samtools
#     - pigz
#     - bcftools
#     - htslib
#     - bedtools
#     - perl
#     - tabix
RUN mkdir -p /conda-envs/4b68e41282e18d46ed0c01f28070afa8
COPY workflow/envs/vanilla_v0.1.yaml /conda-envs/4b68e41282e18d46ed0c01f28070afa8/environment.yaml

# Conda environment:
#   source: workflow/envs/verifybamid2_v0.1.yaml
#   prefix: /conda-envs/a1275c9b142a843967eea08c95835084
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - r

#     - defaults
#   dependencies:
#     - verifybamid2=2.0.1
RUN mkdir -p /conda-envs/a1275c9b142a843967eea08c95835084
COPY workflow/envs/verifybamid2_v0.1.yaml /conda-envs/a1275c9b142a843967eea08c95835084/environment.yaml

# Step 2: Generate conda environments

RUN conda env create --prefix /conda-envs/e332402e4f0181c4533d4808ed2bdab5 --file /conda-envs/e332402e4f0181c4533d4808ed2bdab5/environment.yaml && \
    conda env create --prefix /conda-envs/2c1fcdef6fd8b78abb453cf367d33af6 --file /conda-envs/2c1fcdef6fd8b78abb453cf367d33af6/environment.yaml && \
    conda env create --prefix /conda-envs/1f6a7e34ff1bf45e82e319f82e1f3908 --file /conda-envs/1f6a7e34ff1bf45e82e319f82e1f3908/environment.yaml && \
    conda env create --prefix /conda-envs/d9ae57b0259467e56b5ec3edd572f378 --file /conda-envs/d9ae57b0259467e56b5ec3edd572f378/environment.yaml && \
    conda env create --prefix /conda-envs/13d6c08bdbd8adaab139b9fd6e0ca041 --file /conda-envs/13d6c08bdbd8adaab139b9fd6e0ca041/environment.yaml && \
    conda env create --prefix /conda-envs/eafc04c861b77d27fd56364ce234b91c --file /conda-envs/eafc04c861b77d27fd56364ce234b91c/environment.yaml && \
    conda env create --prefix /conda-envs/4076a60bcc18d220034a09af94d1ce9d --file /conda-envs/4076a60bcc18d220034a09af94d1ce9d/environment.yaml && \
    conda env create --prefix /conda-envs/a4543ea14359b1714e597b0e70c6a64b --file /conda-envs/a4543ea14359b1714e597b0e70c6a64b/environment.yaml && \
    conda env create --prefix /conda-envs/753ea3fe000cd80431db3f4fae0ccbe2 --file /conda-envs/753ea3fe000cd80431db3f4fae0ccbe2/environment.yaml && \
    conda env create --prefix /conda-envs/7b17fcd547d9025990813fe914db3286 --file /conda-envs/7b17fcd547d9025990813fe914db3286/environment.yaml && \
    conda env create --prefix /conda-envs/9bde885de7e411ded58cabb1dae00fb9 --file /conda-envs/9bde885de7e411ded58cabb1dae00fb9/environment.yaml && \
    conda env create --prefix /conda-envs/4e1fe845949f2d4e6edd63dc1c2a3793 --file /conda-envs/4e1fe845949f2d4e6edd63dc1c2a3793/environment.yaml && \
    conda env create --prefix /conda-envs/f080e4234b0c70dc690b6c52fbfe4fbe --file /conda-envs/f080e4234b0c70dc690b6c52fbfe4fbe/environment.yaml && \
    conda env create --prefix /conda-envs/244b237b7cd6c4b6fef25b5fcd603b7a --file /conda-envs/244b237b7cd6c4b6fef25b5fcd603b7a/environment.yaml && \
    conda env create --prefix /conda-envs/e40c66cc926b5d83dd51ff71746f8a18 --file /conda-envs/e40c66cc926b5d83dd51ff71746f8a18/environment.yaml && \
    conda env create --prefix /conda-envs/4b68e41282e18d46ed0c01f28070afa8 --file /conda-envs/4b68e41282e18d46ed0c01f28070afa8/environment.yaml && \
    conda env create --prefix /conda-envs/a1275c9b142a843967eea08c95835084 --file /conda-envs/a1275c9b142a843967eea08c95835084/environment.yaml && \
    conda clean --all -y
