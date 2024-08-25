#!/bin/bash
 
 
 
cp config/MGI_example/procs.txt .
workflow_dir=$PWD
cd $workflow_dir


java -jar /usr/local/bin/womtool.jar validate config/CROMWELL/immuno/analysis-wdls/definitions/immuno.wdl

#java -Dbackend.default=Local -Dconfig.file=config/CROMWELL/app.conf  -jar /usr/local/bin/cromwell.jar run config/CROMWELL/immuno/analysis-wdls/definitions/immuno.wdl -i config/CROMWELL/immuno/example_immuno_cloud-WDL.yaml

#java -Dbackend.default=Local -Dconfig.file=config/MGI_example/app.conf -jar /usr/local/bin/cromwell.jar run config/MGI_example/example.wdl -i config/MGI_example/inputs.json
