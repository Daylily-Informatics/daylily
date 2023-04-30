#!/bin/bash

exec 3>&1


result=$(dialog --inputbox test 0 0 2>&1 1>&3)

echo $result
