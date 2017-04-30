#!/bin/sh

if [ $# != 3 ]; then
  echo "$0 <vcf file> <sample list> <output prefix>"
  exit
fi

bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%TYPE[\t%IUPACGT]\n' -S ${2} ${1} | vcfToSnp > ${3}.snp
