#!/bin/bash
 
 # using the WDLs from https://github.com/wustl-oncology
 
 
cp config/MGI_example/procs.txt .
workflow_dir=$PWD
cd $workflow_dir

## Validate
# java -jar /usr/local/bin/womtool.jar validate config/CROMWELL/immuno/analysis-wdls/definitions/immuno.wdl

## Generate Input JSON Template
# java -jar /usr/local/bin/womtool.jar inputs config/CROMWELL/immuno/analysis-wdls/definitions/immuno.wdl > config/CROMWELL/immuno/example_immuno_cloud-WDL.json    

## Generate Graph
# java -jar /usr/local/bin/womtool.jar graph config/CROMWELL/immuno/analysis-wdls/definitions/immuno.wdl > immuno.dot 
# dot -Tpng immuno.dot -o immuno_wf.png


# Germine WGS
# java -Dbackend.default=local -Dconfig.file=config/CROMWELL/app.conf  -jar /usr/local/bin/cromwell.jar run config/CROMWELL/immuno/analysis-wdls/definitions/germline_wgs.wdl -i config/CROMWELL/immuno/jem_germline_wgs_in.json

# IMMUNO
#java -Dbackend.default=local -Dconfig.file=config/CROMWELL/app.conf  -jar /usr/local/bin/cromwell.jar run config/CROMWELL/immuno/analysis-wdls/definitions/immuno.wdl -i config/CROMWELL/immuno/example_immuno_cloud-WDL.json