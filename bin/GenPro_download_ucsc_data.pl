#!/usr/bin/env perl

# Name: GenPro_download_ucsc_data
#
# Description: Download UCSC data needed to make indexed genome required for
#              SeqAnt and personal protein database work
#
# Date created: Fri Aug 22 14:01:27 2014
# Date last modified: 2016-04-30
# By: TS Wingo

use 5.10.0;
use strict;
use warnings;

use Data::Dump qw/ dump /;
use Getopt::Long;
use Path::Tiny;
use IO::Uncompress::Gunzip qw/ $GunzipError /;

our $VERSION = '0.01';

my @chrs = map { "chr$_" } ( 1 .. 22, "M", "X", "Y" );
my ( $dir_name, $act, $verbose, $genome_name );

die
  "Usage: $0 [--act] [--verbose] --dir <target_dir> --genome <genome to download>\n"
  unless GetOptions(
  'verbose|v'  => \$verbose,
  'act|a'      => \$act,
  'dir|d=s'    => \$dir_name,
  'genome|g=s' => \$genome_name
  )
  and $dir_name
  and $genome_name;

$verbose++ unless $act;

# UCSC places the chromsome field in different positions depending on the
# dataset dump; check the sql dump file for correct order
my @snp_fields = qw/ bin chrom chromStart chromEnd name score strand refNCBI
  refUCSC observed molType class valid avHet avHetSE func locType weight
  exceptions submitterCount submitters alleleFreqCount alleles alleleNs
  alleleFreqs bitfields/;
my @knownGene_fields = qw/ name chrom strand txStart txEnd cdsStart cdsEnd
  exonCount exonStarts exonEnds proteinID alignID/;
my @refGene_fields = qw/bin name chrom strand txStart txEnd cdsStart cdsEnd
  exonCount exonStarts exonEnds score name2 cdsStartStat cdsEndStat exonFrames/;
# NOTE: the kgXref table has 2 more fields - rfamAcc & tRnaName; I don't care
#       about them; and, if those fields are not present then there are no
#       tabs written when UCSC dumps the sql table making them hard to process
#       without just leaving them off the list
my @kgXref_fields =
  qw/ kgID mRNA spID spDisplayID geneSymbol refseq protAcc description /;
my %fields_for_set = (
  'knownGene' => \@knownGene_fields,
  'refGene'   => \@refGene_fields,
  'snp146'    => \@snp_fields,
  'snp142'    => \@snp_fields,
);

# Not all of the UCSC data is needed and we also want to unify the refGene and
# knownGene fields so they can be processed in the same manner
my @wanted_snp_fields = qw/ chrom chromStart chromEnd name alleleFreqCount
  alleles alleleNs alleleFreqs/;
my @wanted_knowGene_fields = qw/ chrom strand txStart txEnd cdsStart cdsEnd
  exonCount exonStarts exonEnds name /;
my @wanted_refGene_fields = qw/ chrom strand txStart txEnd cdsStart cdsEnd
  exonCount exonStarts exonEnds name /;
my @wanted_kgXref_fields =
  qw/ kgID mRNA spID spDisplayID geneSymbol refseq protAcc /;
my %wanted_fields_for_set = (
  'knownGene' => \@wanted_knowGene_fields,
  'refGene'   => \@wanted_refGene_fields,
  'snp146'    => \@wanted_snp_fields,
  'snp142'    => \@wanted_snp_fields,
);

# make base directories
my $db_dir = path($dir_name);
if ( !$db_dir->is_dir && $act ) {
  $db_dir->mkpath;
}

# Rsync flags
my $opt = GetRsyncFlags( $act, $verbose );

# sequence data, per chromosome fasta
dlFromUcsc( $genome_name, $opt, $db_dir->child("chr"), 'chromosomes/chr*.fa.gz' );
if ($act) {
  procFasta( $genome_name, $db_dir->child("chr"), \@chrs );
}
else {
  say "... dry run. use --act to download and process";
}

# cross ref stuff
dlFromUcsc( $genome_name, $opt, $db_dir->child("crossRef"),
  "database/kgXref.txt.gz" );

# gene data, per genome sql dump
#my $gene_set = "knownGene";
my $gene_set = "refGene";
dlFromUcsc( $genome_name, $opt, $db_dir->child("gene"),
  "database/$gene_set.txt.gz" );
if ($act) {
  procGeneSnp( $genome_name, $gene_set, $db_dir->child("gene"),
    \@chrs, $db_dir->child("crossRef/kgXref.txt.gz") );
}
else {
  say "... dry run. use --act to download and process";
}

# snp data, per genome sql dump
#my $snp_set = "snp146";
#dlFromUcsc( $genome_name, $opt, $db_dir->child("snp"), "database/$snp_set.txt.gz" );
#if ($act) {
#  procGeneSnp( $genome_name, $snp_set, $db_dir->child("snp"), \@chrs );
#}
#else {
#  say "... dry run. use --act to download and process";
#}

sub dlFromUcsc {
  my $genome     = shift;
  my $rsyncOpt   = shift;
  my $localDir   = shift;
  my $ucscTarget = shift;

  if ( !$localDir->is_dir ) {
    $localDir->mkpath;
  }
  my $cmd = sprintf( "rsync %s rsync://hgdownload.cse.ucsc.edu/goldenPath/%s/%s %s/",
    $rsyncOpt, $genome, $ucscTarget, $localDir->stringify() );
  say $cmd if $verbose;
  system($cmd) if $act;
}

# procGeneSnp prepares the chromosome-specifc files used by ProcessRefGene.pl
# using the compressed sql dumps available on UCSC's website
sub procGeneSnp {
  my $genome_name    = shift;
  my $set_name       = shift;
  my $dir            = shift;
  my $chrs_aref      = shift;
  my $cross_ref_file = shift;

  # cross reference data for Genes only
  my $cross_ref_href = {};
  if ( defined $cross_ref_file ) {
    $cross_ref_href = procCrossRef( $cross_ref_file, $set_name );
  }

  # the gene and snp files provided by UCSC have different columns for the
  # chromosome; they also do not have any header information; here, we select
  # the right column for the different data sources
  my $field_names_aref = $fields_for_set{$set_name};
  if ( !defined $field_names_aref ) {
    die "unrecognized set '$set_name', could add to 'chr_fields_for_set";
  }
  my $wanted_fields_aref = $wanted_fields_for_set{$set_name};

  # get set name and open the compressed file - all files from UCSC are
  # expected to be gzipped
  my $set_file = $dir->child("$set_name.txt.gz");

  if ($verbose) {
    say "processing '" . $set_file->stringify() . "'...";
  }
  if ( !$act ) {
    if ($verbose) {
      say "done. This is a dry run.";
    }
    return;
  }

  my $z = new IO::Uncompress::Gunzip $set_file->stringify()
    or die "gunzip failed on $set_file: $GunzipError\n";

  # create output files for wanted chromosomes
  my %fh_for_chr =
    map { $_ => $dir->child("$genome_name.$_.$set_name.txt")->filehandle(">") }
    (@$chrs_aref);

  # write data to an allowable chromosome
  my %header = map { $field_names_aref->[$_] => $_ } ( 0 .. $#{$field_names_aref} );
  while (<$z>) {
    chomp $_;
    my @fields = split /\t/, $_;
    if ( scalar keys %header != @fields ) {
      my $msg = sprintf(
        "Error: expected %d fields but found %d",
        scalar keys %header,
        scalar @fields
      );
      die $msg;
    }
    my %data = map { $_ => $fields[ $header{$_} ] } ( keys %header );
    my $extra_info = $cross_ref_href->{ $data{name} };

    if ( exists $fh_for_chr{ $data{chrom} } ) {
      my @out_data = map { $data{$_} } @$wanted_fields_aref;
      if ( defined $extra_info ) {
        push @out_data, $extra_info;
      }
      say { $fh_for_chr{ $data{chrom} } } join( "\t", @out_data );
    }
  }
  if ($verbose) {
    say "done.";
  }
}

sub procCrossRef {
  my $file     = shift;
  my $set_name = shift;

  my ( %info_for_gene, $wanted_col );

  if ( $set_name eq "knownGene" ) {
    $wanted_col = "kgID";
  }
  elsif ( $set_name eq "refGene" ) {
    $wanted_col = "mRNA";
  }
  else {
    die "unknown set: $set_name";
  }
  if ($verbose) {
    my $msg = sprintf( "reading '%s'...", $file->stringify() );
    say $msg;
  }
  my $z = new IO::Uncompress::Gunzip $file->stringify()
    or die "gunzip failed on $file: $GunzipError\n";
  my %header = map { $kgXref_fields[$_] => $_ } ( 0 .. $#kgXref_fields );
  while (<$z>) {
    chomp $_;
    my @fields = split /\t/, $_;
    my %data = map { $_ => $fields[ $header{$_} ] } ( keys %header );
    my @out_data = map { $data{$_} } (@wanted_kgXref_fields);
    $info_for_gene{ $data{$wanted_col} } = join ";", @out_data;
  }
  if ($verbose) {
    say "done";
  }
  return \%info_for_gene;
}

# procFasta uncompresses the chromsome fasta files for wanted chromsomes; this
# is needed for ProcessRefGene.pl
sub procFasta {
  my $genome_name = shift;
  my $localDir    = shift;
  my $chrsAref    = shift;

  for my $chr (@$chrsAref) {
    my $oldName = "$chr.fa.gz";
    my $newName = "$genome_name.$chr.fa";

    my $oldFile = $localDir->child($oldName);
    if ( !$oldFile->is_file && $act ) {
      die "expected file: " . $oldFile->stringify() . " not found\n";
    }
    my $newFile = $localDir->child($newName);
    if ( $newFile->is_file ) {
      next;
    }
    my $cmd =
      sprintf( "gunzip -c %s > %s", $oldFile->stringify(), $newFile->stringify() );
    say $cmd if $verbose;
    system $cmd if $act;
  }
}

sub GetRsyncFlags {
  my $act     = shift;
  my $verbose = shift;

  my $opt;
  # set flags for rsync
  if ($act) {
    if ($verbose) {
      $opt = "-avzP";
    }
    else {
      $opt = "-az";
    }
  }
  else {
    if ($verbose) {
      $opt = "-navzP";
    }
    else {
      $opt = "-naz";
    }
  }
  return $opt;
}
