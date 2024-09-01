## To plug in
## https://github.com/GATB/DiscoSnp


rule discosnp:
    input:
        "nofile.log",
    conda:
        config["discosnp"]["env_yaml"]
    params:
        cluster_sample=ret_sample,
    output:
        "missingfile.log",
    shell:
        "sleep 1;"
