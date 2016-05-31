#!/usr/bin/env perl
# Name:           GenPro_checkPerDb
# Date Created:   Sat May 7 11:33:24 2016
# Date Modified:  Sat May 7 11:33:24 2016
# By:             TS Wingo
#
# Description:    Check dbs created by `make_perprotdb1.pl`

use 5.10.0;
use strict;
use warnings;

use Cpanel::JSON::XS;
use DB_File;
use Data::Dump qw/ dump /;
use Getopt::Long;
use Path::Tiny;

our $VERSION = '0.01';

my ( $id, $db_dir, $verbose );

die "Usage: $0 -d <db dir> -i <id>\n"
  unless GetOptions(
  'v|verbose' => \$verbose,
  'd|dir=s'   => \$db_dir,
  'i|id=s'    => \$id,
  )
  and $id
  and $db_dir;

my $path_db = path($db_dir);
if ( !$path_db->is_dir() ) {
  my $msg = sprintf( "'%s' is not a dir", $db_dir );
  say STDERR $msg;
  exit(1);
}

my @files = path($db_dir)->children(qr/\A$id/);
if ( !@files ) {
  say STDERR "No files found for '$id' in $db_dir";
  exit(1);
}

my ( $gSiteCount, $recCount, $txCount ) = ( 0, 0, 0 );

for my $f (@files) {
  my ( $thisRecCount, $thisTxCount ) = ( 0, 0 );
  my (%chrPos);
  my %db;
  tie %db, 'DB_File', $f->stringify();

  while ( my ( $key, $val ) = each %db ) {
    my $recs_aref = decode_json($val);
    for my $rec_href (@$recs_aref) {
      my $chr = $rec_href->{chr};
      my $pos = $rec_href->{chr_pos};
      $chrPos{"$chr:$pos"}++;
    }
    if ($verbose) {
      say dump( { $key => $recs_aref } );
    }
    $thisRecCount += scalar @$recs_aref;
    $thisTxCount++;
  }
  untie %db;
  printf(
    "File: %s, Replacement sites: %d, transcirpts: %d, record count: %d\n",
    $f, scalar keys %chrPos,
    $thisTxCount, $thisRecCount
  );
  $recCount   += $thisRecCount;
  $txCount    += $thisTxCount;
  $gSiteCount += scalar keys %chrPos;
}

printf( "Total: Replacement sites: %d, transcirpts: %d, record count: %d\n",
  $gSiteCount, $txCount, $recCount );
