#!/usr/bin/env perl
# Name:           GenPro_make_perprotdb1
# Description:    reads the seqant db and a snpfile and prepares a personalized
#                 protein database; use in conjunction with make_refprotdb.pl,
#                 create_db.pl and download_ucsc_data.pl
# Date Created:   Wed Aug 27 21:09:45 2014
# Date Modified:  2017-02-18
# By:             TS Wingo

use 5.10.0;
use strict;
use warnings;

use Cpanel::JSON::XS;
use DB_File;
use Fcntl qw(:DEFAULT :seek);
use Getopt::Long;
use Path::Tiny;
use GenPro qw/ CleanChr /;

our $VERSION = '0.01';

my (
  %db,      %hIUPAC,  $db_name,      $path_bin,   $barcode_file, $snp_file,
  $verbose, $out_dir, $wantedIdFile, $wanted_chr, $debug,
);

#<<< No perltidy
my %codon_2_aa = (
  "AAA" => "K", "AAC" => "N", "AAG" => "K", "AAT" => "N",
  "ACA" => "T", "ACC" => "T", "ACG" => "T", "ACT" => "T",
  "AGA" => "R", "AGC" => "S", "AGG" => "R", "AGT" => "S",
  "ATA" => "I", "ATC" => "I", "ATG" => "M", "ATT" => "I",
  "CAA" => "Q", "CAC" => "H", "CAG" => "Q", "CAT" => "H",
  "CCA" => "P", "CCC" => "P", "CCG" => "P", "CCT" => "P",
  "CGA" => "R", "CGC" => "R", "CGG" => "R", "CGT" => "R",
  "CTA" => "L", "CTC" => "L", "CTG" => "L", "CTT" => "L",
  "GAA" => "E", "GAC" => "D", "GAG" => "E", "GAT" => "D",
  "GCA" => "A", "GCC" => "A", "GCG" => "A", "GCT" => "A",
  "GGA" => "G", "GGC" => "G", "GGG" => "G", "GGT" => "G",
  "GTA" => "V", "GTC" => "V", "GTG" => "V", "GTT" => "V",
  "TAA" => "*", "TAC" => "Y", "TAG" => "*", "TAT" => "Y",
  "TCA" => "S", "TCC" => "S", "TCG" => "S", "TCT" => "S",
  "TGA" => "*", "TGC" => "C", "TGG" => "W", "TGT" => "C",
  "TTA" => "L", "TTC" => "F", "TTG" => "L", "TTT" => "F"
);

$hIUPAC{K}{G} = "T";
$hIUPAC{K}{T} = "G";
$hIUPAC{M}{A} = "C";
$hIUPAC{M}{C} = "A";
$hIUPAC{R}{A} = "G";
$hIUPAC{R}{G} = "A";
$hIUPAC{S}{C} = "G";
$hIUPAC{S}{G} = "C";
$hIUPAC{W}{A} = "T";
$hIUPAC{W}{T} = "A";
$hIUPAC{Y}{C} = "T";
$hIUPAC{Y}{T} = "C";
#>>>

die "Usage: $0
    --s|snpfile <snp_file>
    --l|location <binary genome directory>
    --n|name <genome name (e.g., hg19, mm10, etc.)>
    --o|out <output directory>
  [ --c|chr <chrom name> ]
  [ --b|barcode <barcode Id Lookup (e.g., SL99999 => myId1> ]
  [ --w|wanted <file with Ids to include]\n"
  unless GetOptions(
  's|snpfile=s'  => \$snp_file,
  'b|barcode=s'  => \$barcode_file,
  'n|db_name=s'  => \$db_name,
  'l|location=s' => \$path_bin,
  'v|verbose'    => \$verbose,
  'o|out_dir=s'  => \$out_dir,
  'w|wanted=s'   => \$wantedIdFile,
  'c|chr=s'      => \$wanted_chr,
  )
  and $path_bin
  and $db_name
  and $snp_file
  and $out_dir;

# chr
if ( defined $wanted_chr ) {
  $wanted_chr = CleanChr($wanted_chr);
  Log( "Wanted Chromosome:", $wanted_chr );
}

# check - binary index file location
$path_bin = path($path_bin);
if ( !$path_bin->is_dir() ) {
  Log( "Error", "directory expected for location" );
}

# check - output directory
my $path_out = path($out_dir);
if ( !$path_out->is_dir() ) {
  $path_out->mkpath();
}

# read altnerative ids
my $altIdForId_href = {};
if ( defined $barcode_file ) {
  $altIdForId_href = ReadBarcode($barcode_file);
  Log( "Finished reading", $barcode_file );
}

# read list of ids to print
my $wantedIdHref = {};
if ( defined $wantedIdFile ) {
  $wantedIdHref = ReadWantedIdFile($wantedIdFile);
  Log( "Finished reading", $wantedIdFile );
}

# read snpfile
Log( "Started reading", $snp_file );
my ( $variant_site_href, $ref_for_site_href, $ids_aref ) =
  ReadSnpFile( $snp_file, $wanted_chr, $wantedIdHref, $altIdForId_href );
Log( "Finished reading", $snp_file );

# read binary db to determine replacement sites and write list of sites per chrom
# for each id in snpfile
Log("Started determining replacement sites");
ReplacementSites( $path_bin, $path_out, $variant_site_href, $ref_for_site_href,
  $ids_aref );
Log("Finished determining replacement sites");

# ReplacementSites takes the binary database location and hash references for
# variant information and reference sites, and array reference of ids; it creates
# databases of repacement variants on a per id and per chrom basis in the
sub ReplacementSites {
  my $path_bin     = shift;
  my $path_out     = shift;
  my $variant_href = shift;
  my $ref_href     = shift;
  my $ids_aref     = shift;

  if ( scalar @$ids_aref > 250 ) {
    my $msg = sprintf( "Cannot create personal dbs for >250 and was given '%d'.",
      scalar @$ids_aref );
    Log( "Error", $msg, "Select fewer individuals with --wanted <file with ids>" );
  }

  foreach my $chr ( sort keys %$variant_href ) {

    PrepareDb( $path_out, $chr, $ids_aref );
    Log( "Info", "Prepared Databases", $chr );

    # load the data base, per chrom
    my $fpIDX = $path_bin->child("$db_name.$chr.idx")->filehandle("<");
    binmode $fpIDX;
    my $idx_typedef      = "L S L S C";
    my $idx_typedef_size = length( pack($idx_typedef) );

    my $fpDAT = $path_bin->child("$db_name.$chr.bgd")->filehandle("<");
    binmode $fpDAT;
    my $dat_typedef      = "I C C S C Z25 h40";
    my $dat_typedef_size = length( pack($dat_typedef) );

    my $fpGID = $path_bin->child("$db_name.$chr.gid")->filehandle("<");
    binmode $fpGID;
    my $gid_typedef      = "Z120";
    my $gid_typedef_size = length( pack($gid_typedef) );

    foreach my $pos ( sort { $a <=> $b } keys %{ $variant_href->{$chr} } ) {

      my $min_alleles_aref = $variant_href->{$chr}{$pos}{min_allele};
      my $these_ids_aref   = $variant_href->{$chr}{$pos}{ids};
      my $ref_allele       = $ref_href->{$chr}{$pos};

      my %min_alleles = map { $_ => 1 } (@$min_alleles_aref);
      my @min_alleles = map { $_ } ( sort keys %min_alleles );

      for my $min_allele (@min_alleles) {
        my $idx_offset = ( $pos - 1 ) * $idx_typedef_size;
        seek( $fpIDX, $idx_offset, SEEK_SET );
        if ( read( $fpIDX, my $idx_buffer, $idx_typedef_size ) == $idx_typedef_size
          || die "error reading idx at '$pos'\n" )
        {
          my (
            $dat_file_offset, $dat_record_len, $snp_file_offset,
            $snp_record_len,  $coding_status
          ) = unpack( $idx_typedef, $idx_buffer );

          # if ($verbose) {
          #   my $msg = sprintf(
          #     "dat_offset: %x, dat_len: %x, snp_offset: %x, snp_len: %x, coding_status: %x",
          #     $dat_file_offset, $dat_record_len, $snp_file_offset,
          #     $snp_record_len,  $coding_status
          #   );
          #   say $msg;
          # }

          # skip non-coding sites
          if ( $coding_status != 0 ) {
            next;
          }

          # read dat file
          if ( $dat_record_len > 0 ) {
            my ($new_base);
            seek( $fpDAT, $dat_file_offset, SEEK_SET );
            while ( $dat_record_len > 0 ) {
              my (
                @codon,      @new_codon, $strand,     $gene_symbol, $codon_pos,
                $codon_code, $codon,     $aa_residue, $codon_number,
              );

              read( $fpDAT, my $dat_buffer, $dat_typedef_size ) == $dat_typedef_size
                || die "error reading bgd file at '$dat_file_offset'\n";
              my @gene_dat = unpack( $dat_typedef, $dat_buffer );

              # dat format (here, @gene_dat):
              #      [0] Gene Number
              #      [1] Coding/UTR & Orientation (number 0-9 as key for %genomic_region_codes)
              #      [2] Codon
              #      [3] Amino Acid Code
              #      [4] Error Code
              #      [5] RefSeq Id
              #      [6] sha1 hex -> used b/c RefSeq ID isn't really unique

              # >>> NOT REALLY NEEDED
              # gene name
              my $gid_file_offset = ( $gene_dat[0] - 1 ) * $gid_typedef_size;
              seek( $fpGID, $gid_file_offset, SEEK_SET );
              read( $fpGID, my $gid_buffer, $gid_typedef_size ) == $gid_typedef_size
                || die "error reading gid at $gid_file_offset\n";
              my $gene_desc = unpack( $gid_typedef, $gid_buffer );
              # get individual elements from gene_desc
              my ( $kgID, $mRNA, $spID, $spDisplayID, $geneSymbol, $refseq, $protAcc ) =
                split /\;/, $gene_desc;
              # <<< NOT REALLY NEEDED

              $strand = ( $gene_dat[1] < 5 ) ? '+' : '-';

              # codon and position within codon
              # NOTE: codon is always with respect to the coding strand
              if ( $gene_dat[2] > 0 ) {
                $codon_pos = ( int( ( $gene_dat[2] - 1 ) / 64 ) ) + 1;
                $codon_code = ( $gene_dat[2] - 1 ) % 64;
                for ( my $i = 16; $i >= 1; $i /= 4 ) {
                  $codon .= int( $codon_code / $i );
                  $codon_code = $codon_code % $i;
                }
                $codon =~ tr/0123/ACGT/;
                @codon = split( //, $codon );

                if ( !okCodonRef( $strand, \@codon, $codon_pos, $ref_allele ) ) {
                  my $msg = sprintf( "Wrong reference for %s:%s (+), found '%s' expected '%s'\n",
                    $chr, $pos, $ref_allele, $codon[ $codon_pos - 1 ] );
                  Log( "Info", $msg );
                }
                else {
                  $new_base = $min_allele;
                }
              }

              # AA and AA position
              if ( $gene_dat[3] > 0 ) {
                $aa_residue = $gene_dat[3];
                $codon_number =
                  int( ( $gene_dat[3] * 3 ) + int( ( $gene_dat[2] - 1 ) / 64 ) + 1 - 3 );
              }
              my $err_code        = $gene_dat[4];
              my $transcript_name = $gene_dat[5];
              my $sha1            = $gene_dat[6];

              if ( defined $new_base and defined $aa_residue and defined $strand ) {
                @new_codon = @codon;
                $new_codon[ $codon_pos - 1 ] = $new_base;
                my $new_codon = join "", @new_codon;
                my $new_aa    = $codon_2_aa{$new_codon};
                my $old_aa    = $codon_2_aa{$codon};

                # save the site for certain ids
                if ( defined $new_aa and defined $old_aa and $new_aa ne $old_aa ) {
                  # cycle through the list of IDs and minor alleles and save the transcript and
                  # variant info for each individual with variants in the particular transcripts
                  #                                                   AA_residue    original_AA         new_AA
                  # e.g., push @{ r_sites{$id}{$transcript_name} }, "$aa_residue:$codon2_aa{$codon}:$codon_2_aa{$new_codon}";
                  for ( my $i = 0; $i < @$min_alleles_aref; $i++ ) {
                    my $id        = $these_ids_aref->[$i];
                    my $id_allele = $min_alleles_aref->[$i];
                    if ( $id_allele eq $min_allele ) {
  my $record_href = {
    aa_residue => $aa_residue,
    old_aa     => $codon_2_aa{$codon},
    new_aa     => $codon_2_aa{$new_codon},
    chr        => $chr,
    chr_pos    => $pos,
    codon_pos  => $codon_pos,
    ref_allele => $ref_allele,
    min_allele => $min_allele,
  };

  WriteToDb( $id, $chr, $transcript_name, $record_href );
                    }
                  }
                }
              }
              $dat_record_len -= $dat_typedef_size;
            }
          }
        }
      }
    }
    if ( CloseDb( $chr, $ids_aref ) ) {
      my $msg = sprintf( "done with %s.", $chr );
      Log( "Info", $msg );
    }
    else {
      my $msg = springf( "Error closing %s.", $chr );
      Log( "Error", $msg );
    }
    close $fpIDX;
    close $fpDAT;
    close $fpGID;
  }
}

sub CloseDb {
  my ( $chr, $ids_aref ) = @_;
  if ( !defined $chr || !defined $ids_aref ) {
    return;
  }
  for my $id ( keys %db ) {
    untie %{ $db{$id}{$chr} };
    undef $db{$id}{$chr};
  }
  return 1;
}

sub PrepareDb {
  my ( $path, $chr, $ids_aref ) = @_;
  for my $id (@$ids_aref) {
    if ( defined $id ) {
      my %this_db;
      tie %this_db, 'DB_File', $path->child("$id.$chr.db")->stringify();
      $db{$id}{$chr} = \%this_db;
    }
  }
}

sub WriteToDb {
  my ( $id, $chr, $key, $href ) = @_;
  my $db_href = $db{$id}{$chr};

  if ( !defined $db_href ) {
    my $msg = sprintf( "Error - no db for %s %s", $id, $chr );
    Log( "Fatal", $msg );
  }
  my $val = $db_href->{$key};
  if ( !defined $val ) {
    $db_href->{$key} = encode_json( [$href] );
  }
  else {
    my $recs_aref = decode_json($val);
    push @$recs_aref, $href;
    $db_href->{$key} = encode_json($recs_aref);
  }
}

sub okCodonRef {
  my ( $strand, $codonAref, $codonPos, $refAllele ) = @_;
  my $new_base;

  if ( $strand eq "-" ) {
    $refAllele =~ tr/ACGT/TGCA/;
  }

  if ( $codonAref->[ $codonPos - 1 ] eq $refAllele ) {
    return 1;
  }

  return;
}

# ReadBarcode takes a barcode file, expected to have 2 columns and returns a
# hash of the 1st column to second column
sub ReadBarcode {
  my $file = shift;

  my %idLookup;

  my $fh = path($file)->filehandle("<");
  while (<$fh>) {
    chomp $_;
    my @fields = split /\t/, $_;
    $idLookup{ $fields[0] } = $fields[1];
  }
  return \%idLookup;
}

# ReadSnpFile takes a snpfile and an optionally populated alternative id lookup
# hash and returns a list of all references for sites and a list of all reference
# alleles, a hash reference of all variant sites for ids, and all ids with variants
# Snpfile format:
#   - invariant fields: Fragment  Position  Reference Alleles Allele_Counts Type
#   - variable fields: 'Id1\tId2'
#       followed by 2 columns per ID sequneced (genotype call using IUPAC coding of bases and probability of the call)
# Leave loop with:
#   %variant_sites{ id }{ chr }{ pos } = "ref_allele:minor_base"
#   @ids -> note this has blanks in it, which is intentional
sub ReadSnpFile {
  my $file           = shift;
  my $wanted_chr     = shift;
  my $wanted_id_href = shift;
  my $id_lookup_href = shift;
  my $okChrsAref     = shift;

  my ( %ids, %variant_sites, %variant_sites_for_id, %ref_for_site );

  # set some expectations for the snpfile
  my $snpHeaderFieldEnd = 5;
  #my @expHeader = qw/Fragment Position Reference Alleles Allele_Counts Type/;

  # defined acceptable chromosomes, which is really dependent on the binary db
  if ( !defined $okChrsAref ) {
    my @chrs = map { "chr$_" } ( 1 .. 22, 'X', 'Y', 'M' );
    $okChrsAref = \@chrs;
  }
  my %okChr = map { $_ => 1 } @$okChrsAref;

  my $fh = path($file)->filehandle("<");

  my $headerLine = <$fh>;
  chomp $headerLine;
  my @headerFields = split '\t', $headerLine;
  my %header = map { $headerFields[$_] => $_ } ( 0 .. $snpHeaderFieldEnd );
  my @ids = ProcessIds( \@headerFields, $wanted_id_href, $id_lookup_href );

  while (<$fh>) {
    chomp $_;
    my @fields = split '\t', $_;
    my %data = map { $_ => $fields[ $header{$_} ] } ( keys %header );
    checkSnpData( \%data, $. );
    my $chr  = $data{Fragment};
    my $pos  = $data{Position};
    my $ref  = $data{Reference};
    my $type = $data{Type};
    $ref_for_site{$chr}{$pos} = $ref;

    # we annotate all SNP records that are in acceptable chromosomes
    if ( $type ne "SNP" || !exists $okChr{$chr} ) {
      next;
    }

    # proc genotyping
    for ( my $i = $snpHeaderFieldEnd + 1; $i < @fields; $i += 2 ) {
      my $id = $ids[$i];

      # skip blank ids
      if ( $id eq "" ) {
        next;
      }

      my $geno = $fields[$i];
      my $prob = $fields[ $i + 1 ];
      if ( !defined $geno || !defined $prob ) {
        Log( "Error", "Id: '$id' has no genotype or probability." );
      }
      if ( $geno eq $ref ) {
        next;
      }
      if ( $prob >= 0.95 ) {
        my $min_allele = $hIUPAC{$geno}{ $data{Reference} };
        if ( defined $min_allele ) {
          $ids{ $ids[$i] }++;
          push @{ $variant_sites{$chr}{$pos}{ids} },        $ids[$i];
          push @{ $variant_sites{$chr}{$pos}{min_allele} }, $min_allele;
        }
      }
    }
  }
  # clean up ids for return
  @ids = map { $_ } sort { $a cmp $b } ( keys %ids );
  return ( \%variant_sites, \%ref_for_site, \@ids );
}

sub ProcessIds {
  my ( $fieldsAref, $wantedIdHref, $altNameHref ) = @_;

  my @ids;

  # using the snpfile header line to get ids - format 'id\t\tid' and we split
  # using '\t' in ReadSnpFile
  for ( my $i = 0; $i < @$fieldsAref; $i++ ) {
    my $id = $fieldsAref->[$i];
    if ( !defined $id ) {
      $ids[$i] = "";
      next;
    }
    if (%$altNameHref) {
      if ( exists $altNameHref->{$id} ) {
        $id = $altNameHref->{$id};
      }
    }
    if (%$wantedIdHref) {
      if ( exists $wantedIdHref->{$id} ) {
        $ids[$i] = $id;
      }
      else {
        $ids[$i] = "";
      }
    }
    else {
      $ids[$i] = $id;
    }
  }
  if (wantarray) {
    return @ids;
  }
  elsif ( defined wantarray ) {
    return \@ids;
  }
  else {
    my $msg = "ProcessIds() should be called in list or scalar context.";
    Log( "Error", $msg );
  }
}

sub checkSnpData {
  my ( $href, $lineNum ) = @_;

  for my $reqField (qw/ Fragment Position Reference Type /) {
    if ( !exists $href->{$reqField} ) {
      my $msg =
        sprintf( "Snpfile missing required field: '%s' at line %d", $reqField, $lineNum );
      Log( "Error", $msg );
    }
    elsif ( !defined $href->{$reqField} ) {
      my $msg =
        sprintf( "Snpfile field: '%s' undefined  at line %d", $reqField, $lineNum );
      Log( "Error", $msg );
    }
  }
}

# ReadWantedIdFile takes a file name and returns a hash reference of a list of
# Ids within the 1st column of the file (tab separated)
sub ReadWantedIdFile {
  my $file = shift;

  my %wantedIds;

  my @lines = path($file)->lines( { chomp => 1 } );
  for my $line (@lines) {
    ( my $clean_line = $line ) =~ s/\s+\z//xm;
    if ( defined $clean_line ) {
      $wantedIds{$clean_line} = 1;
    }
    else {
      my $msg = sprintf( "Error: got '%s' id for line '%s' from file: %s",
        $clean_line, $line, $file );
      Log( "Error", $msg );
    }
  }
  return \%wantedIds;
}

# Log takes and array of strings and uses the first element to handle the
# message reporting - 1) Error, probably a user/input error (print and exit),
# 2) Fatal - probably an internal error, print and croak, 3) Debug is
# something that is printed if the global debug is true, and 4) Warn prints
# a message to standard out regardless of the level of verbosity. Anything else
# will be printed depending on whether the gloabal verbose value is set.
sub Log {
  my $type = shift;
  if ( $type eq 'Error' ) {
    my $msg = join ": ", $type, ( join " ", @_ );
    say STDERR $msg;
    exit(1);
  }
  elsif ( $type eq 'Fatal' ) {
    my $msg = join ": ", $type, ( join " ", @_ );
    croak $msg;
  }
  elsif ( $type eq 'Warn' ) {
    my $msg = join ": ", $type, ( join " ", @_ );
    say STDERR $msg;
  }
  elsif ( $type eq 'Info' ) {
    if ($verbose) {
      my $msg = join ": ", $type, ( join " ", @_ );
      say STDERR $msg;
    }
  }
  elsif ( $type eq 'Debug' ) {
    if ($debug) {
      my $msg = join ": ", $type, ( join " ", @_ );
      say STDERR $msg;
    }
  }
  else {
    if ($verbose) {
      my $msg = join " ", $type, @_;
      say STDERR "Info: ", $msg;
    }
  }
}
