# ### QUIT Imputation and HLA Typing
#
# https://github.com/rwdavies/QUILT#paragraph-installation-conda
# https://www.nature.com/articles/s41588-021-00877-0


rule quilt:
    conda:
        "../envs/quilt_v0.1.yaml"
    output:
        "logs/quilt.done",
    params:
        cluster_sample=f"wildcards.sample",
    shell:
        "touch {output}"
