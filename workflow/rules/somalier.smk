# ##### SOMALIER - A burlier peddy
# --------------------------------
# Infer relatedness, ethinicity, QC applications
# github: https://github.com/brentp/somalier/releases/tag/v0.2.13
# https://hub.docker.com/r/brentp/somalier/


rule somalier:
    params:
        cluster_sample=ret_sample,
    container:
        "docker://brentp/somalier"
    shell:
        "echo  but have fetched and stored all of the source data/binaries"
