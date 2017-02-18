#!/bin/sh

if [ $# != 2 ]; then
  echo "$0 <vcf file> <sample list>"
  exit
fi

bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%TYPE[\t%IUPACGT]\n' -S ${2} ${1} | vcfToSnp > ${1}.snp
