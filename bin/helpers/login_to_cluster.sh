#!/bin/bash
pcluster describe-cluster -n $cluster_name | grep 'publicIpAddress' | perl -pe 's/.*\: \"(.*)\"\,.*/$1/g;' &&  cluster_ip=$(pcluster describe-cluster -n $cluster_name | grep 'publicIpAddress' | perl -pe 's/.*\: \"(.*)\"\,.*/$1/g;') && ssh -i ~/.ssh/daylily2.pem ubuntu@$cluster_ip
