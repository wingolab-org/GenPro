#!/bin/sh
GenPro_create_db.pl --g ${1} --c ${SGE_TASK_ID} --genedir ${1}/gene --geneset knownGene -f ${1}/chr -o ${1}/idx -v &> ${1}.${SGE_TASK_ID}.idxBuild.log
