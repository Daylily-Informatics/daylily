#!/bin/bash

echo $(curl -s http://169.254.169.254/latest/meta-data/instance-type)
