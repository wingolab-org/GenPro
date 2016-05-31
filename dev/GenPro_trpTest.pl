#!/usr/bin/perl

use 5.10.0;
use warnings;
use strict;

use Data::Dump qw/dump/;

my $min_peptide_length = 6;
my $max_peptide_length = 45;
my $trim_end_regex     = qr{\*[\*\w]*\z};

my $prot1 =
  "MWWKQLVAGAVAGAVSRTGTAPLDRLKVFMQVHASKTNRLNILGGLRSMVLEGGIRSLWRGNGINVLKIAPESAIKFMAYEQIKRAILGQQETLHVQERFVAGSLAGATAQTIIYPMETLKNWWLQQYSHDSADPGILVLLACGTISSTCGQIASYPLALVRTRMQAQASIEGGPQLSMLGLLRHILSQEGMRGLYRGIAPNFMKVIPAVSISYVVYENMKQALGVTSR";
my $prot2 =
  "MWWKQLVAGAVAGAVSRTGTAPLDRLKVFMQVHASKTNRLNILGGLRSMVLEGGIRSLWRGNGINVLKIAPESAIKFMAYEQIKRAILGQQETLHVQERFVAGSLAGATAQTIIYPMETLKNWWLQQYSHDSADPGILVLLACGTISSTCGQIASYPLALVRTRMQAQASIEGGPQLSMLGLLRHILSQEGMRGLYRGIAPNFMKVIPAVSISYVVYENMKQALGVTSR";

my %test = (
  ''                     => { Fragments => [], },
  AAAAVVVRRRTTTTTTTTKRRT => { Fragments => [ "AAAAVVVR", "TTTTTTTTK" ], },
  AAAAARCCCCCRP          => { Fragments => [ "AAAAAR", "CCCCCR", "CCCCCRP" ], },
  AAAARPVVVR             => { Fragments => ["AAAARPVVVR"], },
  AAAARPVVVVR => { Fragments => [ "AAAARPVVVVR", "PVVVVR" ], },
  AAAARPVVVVR => { Fragments => [ "AAAARPVVVVR", "PVVVVR" ], },
  $prot1      => {
    Fragments => [
      "QLVAGAVAGAVSR",        "TGTAPLDR",
      "VFMQVHASK",            "LNILGGLR",
      "SMVLEGGIR",            "GNGINVLK",
      "IAPESAIK",             "FMAYEQIK",
      "AILGQQETLHVQER",       "FVAGSLAGATAQTIIYPMETLK",
      "MQAQASIEGGPQLSMLGLLR", "HILSQEGMR",
      "GIAPNFMK",             "VIPAVSISYVVYENMK",
      "QALGVTSR"
    ],
  },
  $prot2 => {
    Fragments => [
      "QLVAGAVAGAVSR",                             "TGTAPLDR",
      "VFMQVHASK",                                 "LNILGGLR",
      "SMVLEGGIR",                                 "GNGINVLK",
      "IAPESAIK",                                  "FMAYEQIK",
      "AILGQQETLHVQER",                            "FVAGSLAGATAQTIIYPMETLK",
      "NWWLQQYSHDSADPGILVLLACGTISSTCGQIASYPLALVR", "MQAQASIEGGPQLSMLGLLR",
      "HILSQEGMR",                                 "GIAPNFMK",
      "VIPAVSISYVVYENMK",                          "QALGVTSR"
    ],
  },
);

for my $seq ( sort keys %test ) {
  my %d = Trypsin($seq);
  say dump( { seq => $seq, digest => \%d, expected_dig => $test{$seq} } );
  checkResults( $test{$seq}{Fragments}, [ values %d ] );
}

sub checkResults {
  my ( $expAref, $obsAref ) = @_;

  my @sExp = sort { $a cmp $b } @$expAref;
  my @sObs = sort { $a cmp $b } @$obsAref;

  my $msg = checkSortedList( \@sExp, \@sObs );
  if ( $msg ne "" ) {
    say STDERR $msg, "\n", dump( { obs => \@sObs, exp => \@sExp } );
    exit(1);
  }
  else {
    say "--- Ok ---";
  }
}

sub checkSortedList {
  my ( $aref1, $aref2 ) = @_;

  my $msg = "";

  if ( @$aref1 != @$aref2 ) {
    $msg = "not equal lengths\n";
  }

  for ( my $i = 0; $i < @$aref1; $i++ ) {
    if ( !exists $aref2->[$i] ) {
      return $msg;
    }
    if ( $aref1->[$i] ne $aref2->[$i] ) {
      $msg .=
        sprintf( "arrays differ at %d element: %s, %s", $i, $aref1->[$i], $aref2->[$i] );
      return $msg;
    }
  }
  return $msg;
}

sub Trypsin {
  my $peptide = shift;

  my %digest;

  # always trim off anything after a '*'; the if ( index ...) > -1) is (probably)
  # speeding things up since most of the time there's no '*' and using the
  # regex is comparitively slow versus just checking whether the string exists
  if ( index( $peptide, '*' ) > -1 ) {
    $peptide =~ s/$trim_end_regex//xm;
  }
  my @peptide = split( //, $peptide );
  my @cut_sites;
  my $last_cut_site = 0;
  for ( my $i = 0; $i < @peptide; $i++ ) {
    if ( $peptide[$i] eq "K" || $peptide[$i] eq "R" ) {
      push @cut_sites, $i + 1;
      $last_cut_site = $i + 1;
    }
  }

  # if a K or R at the end, don't include it twice but also include the
  # end of the protein
  if ( $last_cut_site < length $peptide ) {
    push @cut_sites, scalar @peptide;
  }
  say dump( { cut_sites => \@cut_sites } );

  # reset
  $last_cut_site = 0;

  for ( my $i = 0; $i < @cut_sites; $i++ ) {
    my $start = $last_cut_site;
    my $end   = $cut_sites[$i];
    my $seq   = join "", @peptide[ $start .. $end - 1 ];

    if ( length $seq >= $min_peptide_length && length $seq <= $max_peptide_length ) {
      my $coord = sprintf( "%d:%d", $start + 1, $end );
      $digest{$coord} = $seq;
    }

    if ( @peptide == $end ) {
      last;
    }

    if ( $peptide[$end] eq "P" ) {

      # any more cut sites?
      if ( $i + 1 <= @cut_sites ) {
        my $end = $cut_sites[ $i + 1 ];
        my $seq = join "", @peptide[ $start .. $end - 1 ];
        if ( length $seq >= $min_peptide_length && length $seq <= $max_peptide_length ) {
          my $coord = sprintf( "%d:%d", $start + 1, $end );
          $digest{$coord} = $seq;
        }
      }
    }
    $last_cut_site = $end;
  }
  if (wantarray) {
    return %digest;
  }
  elsif ( defined wantarray ) {
    return \%digest;
  }
  else {
    die "Trypsin() should be called in the list or scalar context";
  }
}
