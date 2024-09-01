import os
import sys

##### MANAGING SUPPORTING DATA
# ----------------------------


localrules: remove_supporting_data, force_supporting_data_cleanup,
rule remove_supporting_data:
    input:
        "logs/all_cleanup.ready",
    output:
        touch("logs/supporting_data_cleanup.done"),
    threads: 1
    log:
        "logs/supporting_data_cleanup.log",
    conda:
        config["vanilla"]["env_yaml"]
    script:
        "../scripts/cleanup_supporting_data.py"


rule force_supporting_data_cleanup:  # TARGET : -not implemented yet- Will Clean Up The Supporing Data Situation If Auto Cleanup Has Failed
    output:
        touch("logs/forced_supporting_data_cleanup.done"),
    log:
        "logs/forced_supporting_data_cleanup.log",
    threads: 1
    conda:
        config["vanilla"]["env_yaml"]
    script:
        "../scripts/cleanup_supporting_data.py"


localrules:
    stage_supporting_data,


rule stage_supporting_data:  # TARGET : Create resources/fsx Link via Method Chosen In Config To Supporting Data DAY Requires (references, etc).
    """ The Action Is Really Happening In the Python Script """
    output:
        "logs/supporting_data_staging.done",
    threads: 1
    log:
        "logs/staging_supporting_data.log",
    conda:
        config["vanilla"]["env_yaml"]
    shell:
        "touch {output};"
        "{latency_wait};"
