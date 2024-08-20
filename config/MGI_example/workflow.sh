#!/bin/bash
 
set -e
err_report() {
  echo "Error on line $1, error code:$?"
}
trap 'err_report $LINENO' ERR
 
 
workflow_dir=/path/to/test-project
cd $workflow_dir


java -Dbackend.default=Local -Dconfig.file=config/MGI_example/app.conf -jar /usr/local/bin/cromwell.jar run config/MGI_example/example.wdl -i config/MGI_example/inputs.json
