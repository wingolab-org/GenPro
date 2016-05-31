#!/bin/sh
./GenPro_create_db.pl --g hg38 --c ${SGE_TASK_ID} --genedir hg38/gene --geneset knownGene -f hg38/chr -o hg38/idx -v &> hg38.chr${SGE_TASK_ID}.idxBuild.log
