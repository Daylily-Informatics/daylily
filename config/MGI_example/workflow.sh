#!/bin/bash
 
 
 
cp config/MGI_example/procs.txt .
workflow_dir=$PWD
cd $workflow_dir


java -Dbackend.default=Local -Dconfig.file=config/MGI_example/app.conf -jar /usr/local/bin/cromwell.jar run config/MGI_example/example.wdl -i config/MGI_example/inputs.json
