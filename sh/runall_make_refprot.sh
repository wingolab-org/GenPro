#!/bin/bash

# $1 -> genome build
# $2 => binary genome dir
# $3 => output dir

./GenPro_make_refprotdb.pl -l ${2} -n ${1} -c chr${SGE_TASK_ID} -o ${3} &> ${1}.chr${SGE_TASK_ID}.build.log
