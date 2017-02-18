#!/usr/bin/env perl
# Name:           GenPro_make_refprotdb
# Description:    reads the seqant db and makes a reference protein database
#                 based on all available transcripts; use in conjuction with
#                 create_db.pl and download_ucsc_data.pl
#
# Date Created:   Sat Aug 23 20:12:01 2014
# Date Modified:  2017-02-18
# By:             TS Wingo

use 5.10.0;
use strict;
use warnings;

use Cpanel::JSON::XS;
use Fcntl qw(:DEFAULT :seek);
use Getopt::Long;
use Path::Tiny;
use GenPro qw/ CleanChr /;

our $VERSION = '0.01';

my ( $chr, $out_path, $idx_path, $db_name, $description_file, $verbose, $debug );
my ( %proteins, %transcripts, %aa_pos, %gene_symb_descipt );

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
my %hCoding_status = ( "0" => "Coding", "1" => "Intergenic", "2" => "Intronic" );
#>>>

die "Usage: $0
    --l|location <binary genome location>
    --n|name <name of seqant db, e.g., hg19, hg18, mm10, etc.>
    --c|chr <chromosome name>
    --o|out <out dir>
  [ --v|verbose ]\n"
  unless GetOptions(
  'l|location=s' => \$idx_path,
  'n|db_name=s'  => \$db_name,
  'o|out=s'      => \$out_path,
  'c|chr=s'      => \$chr,
  'v|verbose'    => \$verbose,
  )
  and $idx_path,
  and $db_name
  and $out_path
  and $chr;

# sanity check
$chr = CleanChr($chr);

$out_path = path($out_path);
if ( !$out_path->is_dir() ) {
  $out_path->mkpath();
}
$idx_path = path($idx_path);
if ( !$idx_path->is_dir() ) {
  die "expected location to be a directory with binary db";
}

# main
say "reading genome for '$chr'...";
ReadDb();
say "done.";

say "Writing protein data for '$chr'...";
WriteRefProt($out_path);
say "done";

# functions

# ReadDb process binary index files to extract all reference proteins
sub ReadDb {
  my $fpIDX = $idx_path->child("$db_name.$chr.idx")->filehandle("<");
  binmode $fpIDX;
  my $idx_typedef      = "L S L S C";
  my $idx_typedef_size = length( pack($idx_typedef) );

  my $fpDAT = $idx_path->child("$db_name.$chr.bgd")->filehandle("<");
  binmode $fpDAT;
  my $dat_typedef      = "I C C S C Z25 h40";
  my $dat_typedef_size = length( pack($dat_typedef) );

  my $fpGID = $idx_path->child("$db_name.$chr.gid")->filehandle("<");
  binmode $fpGID;
  my $gid_typedef      = "Z120";
  my $gid_typedef_size = length( pack($gid_typedef) );

  my $pos = 0;

  while ( read( $fpIDX, my $idx_buffer, $idx_typedef_size ) ) {
    # read index file
    $pos++;
    my (
      $dat_file_offset, $dat_record_len, $snp_file_offset,
      $snp_record_len,  $coding_status
    ) = unpack( $idx_typedef, $idx_buffer );

    # skip non-coding regions
    if ( $coding_status != 0 ) {
      next;
    }

    if ($verbose) {
      my $msg = sprintf(
        "dat_offset: %x, dat_len: %x, snp_offset: %x, snp_len: %x, coding_status: %x",
        $dat_file_offset, $dat_record_len, $snp_file_offset,
        $snp_record_len,  $coding_status
      );
      say "$chr:$pos, $hCoding_status{$coding_status} => $msg";
    }

    # read dat file
    if ( $dat_record_len > 0 ) {
      my $dat_buffer;
      seek( $fpDAT, $dat_file_offset, SEEK_SET );
      while ( $dat_record_len > 0 ) {
        read( $fpDAT, $dat_buffer, $dat_typedef_size ) == $dat_typedef_size
          || die "error reading bgd file at address $dat_file_offset\n";
        my @gene_dat = unpack( $dat_typedef, $dat_buffer );
        my $codon;

        # Gene Dat (i.e., @gene_dat) format:
        #      [0] Gene Number -> offset for gid file (gene symbol)
        #      [1] Coding status and strand  -> 1 - 9 ( FWD strand < 5) coding status
        #      [2] Codon
        #      [3] Amino Acid Position Code
        #      [4] Error Code -> warn if not zero
        #      [5] RefSeq Id -> keys
        #      [6] sha1_hex -> add to refseq ID

        # read GID (gene name)
        my $gid_file_offset = ( $gene_dat[0] - 1 ) * $gid_typedef_size;
        seek( $fpGID, $gid_file_offset, SEEK_SET );
        read( $fpGID, my $gid_buffer, $gid_typedef_size ) == $gid_typedef_size
          || die "error reading gid file address $gid_file_offset\n";
        my $gene_desc = unpack( $gid_typedef, $gid_buffer );

        # strand
        my $strand = ( $gene_dat[1] < 5 ) ? '+' : '-';

        # codon and position within codon
        if ( $gene_dat[2] > 0 ) {
          my $codon_pos = ( int( ( $gene_dat[2] - 1 ) / 64 ) ) + 1;
          my $codon_code = ( $gene_dat[2] - 1 ) % 64;
          for ( my $i = 16; $i >= 1; $i /= 4 ) {
            $codon .= int( $codon_code / $i );
            $codon_code = $codon_code % $i;
          }
          $codon =~ tr/0123/ACGT/;
        }
        my $aa_pos  = $gene_dat[3];
        my $err     = $gene_dat[4];
        my $gene_id = $gene_dat[5];
        my $sha1    = $gene_dat[6];

        # the following is to avoid duplicateing amino acids for a particular NM_transcript
        # because we're going by bp position (NOT codon position)
        if ( defined $codon and $aa_pos > 0 ) {
          die "error with codon: $codon" unless ( exists( $codon_2_aa{$codon} ) );
          unless ( exists( $aa_pos{"$gene_id:$sha1:$strand:$gene_desc"}{$aa_pos} ) ) {
            push @{ $proteins{"$gene_id:$sha1:$strand:$gene_desc"} },    $codon_2_aa{$codon};
            push @{ $transcripts{"$gene_id:$sha1:$strand:$gene_desc"} }, $codon;
          }
          $aa_pos{"$gene_id:$sha1:$strand:$gene_desc"}{$aa_pos}++;
        }
        $dat_record_len -= $dat_typedef_size;
      }
    }
  }
}

# digest proteins, keep peptides that pass size requirements and
# write the hash as a file in JSON format
sub WriteRefProt {
  my $out_path = shift;
  my $json_fh  = path($out_path)->child("refProt.$chr.json")->filehandle(">");
  my $log_fh   = path($out_path)->child("refProt.$chr.log")->filehandle(">");
  my $fasta_fh = path($out_path)->child("refProt.$chr.fasta")->filehandle(">");

  my @tx;
  my ( $good_stop, $bad_stop )    = ( 0, 0 );
  my ( $m_start,   $non_m_start ) = ( 0, 0 );
  foreach my $transcript ( sort keys %proteins ) {
    my ( $gene_id, $sha1, $strand, $gene_desc ) = split( /:/, $transcript );
    my @prot_seq  = @{ $proteins{$transcript} };
    my $prot_seq  = join( "", @prot_seq ) if @prot_seq;
    my @trans_seq = @{ $transcripts{$transcript} };
    my $trans_seq = join( "", @trans_seq ) if @trans_seq;
    my ( $err_start, $err_stop ) = ( 0, 0 );

    # handle transcripts coming from the reverse strand
    if ( $strand eq "-" ) {
      $prot_seq = reverse $prot_seq;
      @prot_seq = split( //, $prot_seq );

      # flipping the transcript is *trickier* than it would seem because it's built from
      # the codons, which are now in reverse order; so, we need to start from the last position
      # and step back every 3 bases and then add the intervening 3 bases to a new array.
      #
      # e.g., Codons (periods placed to make it easier to see)
      #        4   3   2   1
      #       AAT.GCT.GCA.ATG
      my @rev_trans_seq = split( //, $trans_seq );
      my @fwd_trans_seq = ();
      my $codon_number  = ( scalar(@rev_trans_seq) ) / 3;
      for ( my $i = $codon_number; $i > 0; $i-- ) {
        for ( my $j = ( 3 * ( $i - 1 ) ); $j < ( 3 * $i ); $j++ ) {
          push @fwd_trans_seq, $rev_trans_seq[$j];
        }
      }
      @trans_seq = @fwd_trans_seq;
      $trans_seq = join( "", @trans_seq );
    }

    # Does the transcript start with an M and end with a stop ('*')
    ( $prot_seq[0] ne "M" )  ? $err_start++ : undef;
    ( $prot_seq[-1] ne "*" ) ? $err_stop++  : undef;

    #say join("\t", $gene_id, $strand, $gene_symbol, $seq);
    my $err_code = $err_start + $err_stop;
    ($err_start) ? $non_m_start += $err_start : $m_start++;
    ($err_stop)  ? $bad_stop    += $err_stop  : $good_stop++;

    # get individual elements from gene_desc
    my ( $kgID, $mRNA, $spID, $spDisplayID, $geneSymbol, $refseq, $protAcc ) =
      split /\;/, $gene_desc;

    # save relevant data to hash
    my $href = {
      ref_prot    => 1,
      id          => $gene_id,
      desc        => $gene_desc,
      prot_seq    => $prot_seq,
      tx_seq      => $trans_seq,
      strand      => $strand,
      error_code  => $err_code,
      kgID        => $kgID,
      mRNA        => $mRNA,
      spID        => $spID,
      spDisplayID => $spDisplayID,
      geneSymbol  => $geneSymbol,
      refseq      => $refseq,
      protAcc     => $protAcc,
      sha1        => $sha1,
    };
    push @tx, $href;
  }

  # print json
  print $json_fh Cpanel::JSON::XS->new->utf8->encode( \@tx );
  close $json_fh;

  # print log
  say $log_fh
    sprintf( "%d transcripts starting with M, %d transcripts do not start with M",
    $m_start, $non_m_start );
  say $log_fh
    sprintf( "%d transcripts end with stop, %d transcripts do not end with a stop",
    $good_stop, $bad_stop );
  close $log_fh;

  #print fasta
  for my $href (@tx) {
    my $header = sprintf( ">%s |ref_prot;%s", $href->{id}, $href->{desc} );
    print $fasta_fh join( "\n", $header, $href->{prot_seq} ) . "\n";
  }
  close $fasta_fh;
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
