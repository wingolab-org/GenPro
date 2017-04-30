#!/usr/bin/env perl

# NOTE: This is a quick adaptation of Amol Shetty's Old SeqAnt code to created
#       a binary index of the genome useful for creating personal protein dbs.
#       I made some fixes to the original code, mostly in how it handles the
#       construction of the transcript and sequencing sequence. Otherwise, I
#       only made slight modifications, mostly putting things into functions,
#       how the data is packed, and how it expects the gene data, which has
#       been standardized to enable use with knownGene and refGene tracks from
#       UCSC. No further development is expected on this version. -TS Wingo

use 5.10.0;
use strict;
use warnings;

use Carp;
use Digest::SHA qw(sha1_hex);
use Getopt::Long;
use Pod::Usage;
use Path::Tiny;

our $VERSION = '0.01';

##############################################################################
### Constants
##############################################################################

use constant FALSE => 0;
use constant TRUE  => 1;

use constant VERSION => '1.0.0';
use constant PROGRAM => eval { ( $0 =~ m/(\w+\.pl)$/ ) ? $1 : $0 };

use constant FWD_5P_UTR => 0;
use constant FWD_SP_APR => 1;
use constant FWD_CODING => 2;
use constant FWD_SP_DNR => 3;
use constant FWD_3P_UTR => 4;
use constant REV_3P_UTR => 5;
use constant REV_SP_DNR => 6;
use constant REV_CODING => 7;
use constant REV_SP_APR => 8;
use constant REV_5P_UTR => 9;

use constant NON_CODING => 0;
use constant SPLICE_LEN => 6;

use constant PROPER_EXN => 0;
use constant IMPRPR_FRM => 1;
use constant INCMPL_LEN => 2;
use constant BOTH_ERROR => 3;

use constant CODING_POS => 0;
use constant INTERGENIC => 1;
use constant INTRONIC   => 2;

##############################################################################
### Globals
##############################################################################

my %hCmdLineOption = ();
my $sHelpHeader    = "\nThis is " . PROGRAM . " version " . VERSION . "\n";

my ( %hData, %hGeneName, %hdbSnp, %hGeneCoords );
my (@aGenes);
my (
  $sGenome,  $sChr,      $sOutDir,  $sFastaFile, $sInFile,
  $sOutFile, $sGeneFile, $sSnpFile, $sIdxFile
);
my ( $fpREF, $fpFAS, $fpOUT, $fpGID, $fpSNP, $fpIDX );
my ( $sFaHeader, $sFaSequence, $nFaSeqLen );
my (
  $sRefSeqName, $sChrom, $sStrand, $nTStart, $nTEnd,
  $nCStart,     $nCEnd,  $nEStart, $nEEnd
);
my ( $nExonCnt, @aExonStarts,  @aExonEnds,  $sGeneName );
my ( $i,        $n,            $nPos,       $nSplPos );
my ( $sSubSeq,  $nSubSeqStart, $nSubSeqEnd, $nSubSeqLen );
my ( $bDebug,   $bVerbose,     $nError,     $bPreFlag, $bPostFlag );

##############################################################################
### Main
##############################################################################

GetOptions(
  \%hCmdLineOption, 'genedir=s',  'fastadir|f=s', 'dbsnpdir|d=s',
  'geneset=s',      'outdir|o=s', 'genome|g=s',   'chr|c=s',
  'chrstart|s=n',   'chrend|e=n', 'snpbld|n=s',   'bzip2|bz',
  'verbose|v',      'debug',      'help',         'man'
) or pod2usage(2);

if ( $hCmdLineOption{'help'}
  || ( !defined $hCmdLineOption{'geneset'} )
  || ( !defined $hCmdLineOption{'genedir'} )
  || ( !defined $hCmdLineOption{'fastadir'} )
  || ( !defined $hCmdLineOption{'outdir'} )
  || ( !defined $hCmdLineOption{'genome'} )
  || ( !defined $hCmdLineOption{'chr'} ) )
{
  pod2usage( -msg => $sHelpHeader, -exitval => 1 );
}
pod2usage( -exitval => 0, -verbose => 2 ) if $hCmdLineOption{'man'};

$bDebug   = ( defined $hCmdLineOption{'debug'} )   ? TRUE : FALSE;
$bVerbose = ( defined $hCmdLineOption{'verbose'} ) ? TRUE : FALSE;

$sGenome = $hCmdLineOption{'genome'};
$sChr    = CleanChr( $hCmdLineOption{'chr'} );

my $outDir = path( $hCmdLineOption{'outdir'} );
if ( !$outDir->is_dir ) {
  $outDir->mkpath;
}

# prepare files
$sOutFile   = $outDir->child("$sGenome.$sChr.bgd");
$sGeneFile  = $outDir->child("$sGenome.$sChr.gid");
$sSnpFile   = $outDir->child("$sGenome.$sChr.snp");
$sIdxFile   = $outDir->child("$sGenome.$sChr.idx");
$sFastaFile = path( $hCmdLineOption{fastadir} )->child("$sGenome.$sChr.fa");

# fasta
print STDERR "reading fasta file...";
$fpFAS = $sFastaFile->filehandle("<");
ReadFasta($fpFAS);
say STDERR "done.";
$nFaSeqLen = length($sFaSequence);

# Gene
print STDERR "reading gene file...";
my $gene_set = $hCmdLineOption{geneset};
$sInFile = path( $hCmdLineOption{'genedir'} )->child("$sGenome.$sChr.$gene_set.txt");
$fpREF   = $sInFile->filehandle("<");
ReadGene($fpREF);
say STDERR "done.";

# Snp
if ( ( defined $hCmdLineOption{'dbsnpdir'} )
  && ( defined $hCmdLineOption{'snpbld'} ) )
{
  Parse_dbSnp( \%hCmdLineOption, \%hdbSnp );
}
else {
  ( $bDebug || $bVerbose )
    ? print STDERR "Skipping dbSnp data generation .....\n"
    : undef;
}

# wrte binary data
print STDERR "writing binary files...";
BuildBinary();
say "done.";

##############################################################################
### Subroutines
##############################################################################

sub ReadFasta {
  my $fh = shift;
  while (<$fh>) {
    $_ =~ s/\s+$//; # extended chomp

    if ( $_ =~ /^>(.+)/ ) {
      $sFaHeader   = $1;
      $sFaSequence = "";
      next;
    }
    $sFaSequence .= $_;
  }
  close($fh);
}

sub ReadGene {
  my $fh = shift;

  # to allow knownGene and refGene to be used, the download_ucsc_data.pl script
  # re-orders the data, and it tacks on some additional information from kgXref
  # in a standard format: kgID mRNA spID spDisplayID geneSymbol refseq protAcc
  my @expected_fields = qw/ chrom strand txStart txEmd cdsStart cdsEnd
    exonCount exonStarts exonEnds name info/;

  my %header = map { $expected_fields[$_] => $_ } ( 0 .. $#expected_fields );
  my $nX     = 0;
  my $nN     = 0;
  my $nNumAA = 0;

  while (<$fh>) {
    chomp $_;
    my @aFields = split( /\t/, $_ );
    my %data = map { $_ => $aFields[ $header{$_} ] } ( keys %header );

    my ( @this_transcript_seq, @this_transcript_pos, @this_coding_seq, );

    my $sha1_digest = sha1_hex( join( ",", @aFields ) );

    $sChrom      = $data{chrom};
    $sStrand     = $data{strand};
    $nTStart     = $data{txStart};
    $nTEnd       = $data{txEnd};
    $nCStart     = $data{cdsStart};
    $nCEnd       = $data{cdsEnd};
    $nExonCnt    = $data{exonCount};
    @aExonStarts = split( /,/, $data{exonStarts} );
    @aExonEnds   = split( /,/, $data{exonEnds} );
    $sGeneName   = $data{info} || $data{name};
    $sRefSeqName = $data{name};
    $nError      = 0;

    # check validity of RefGene Entry
    if ( scalar @aExonEnds != $nExonCnt || scalar @aExonStarts != $nExonCnt ) {
      my $msg = "Error: count of exons and number of start/stops disagree";
      $msg .= "\n\tLine: $_";
      croak $msg;
    }

    if ( !defined $hGeneName{$sGeneName} ) {
      $nN++;
      $hGeneName{$sGeneName}[0] = $nN;             # Gene Number
      $hGeneName{$sGeneName}[1] = 0;               # Gene Occurrence
      $hGeneName{$sGeneName}[2] = $aFields[4] + 1; # Gene Start
      $hGeneName{$sGeneName}[3] = $aFields[5];     # Gene End
      push @aGenes, $sGeneName;
    }
    $hGeneName{$sGeneName}[1] = $nX;
    $hGeneName{$sGeneName}[2] = ( $aFields[4] + 1 )
      if ( ( $aFields[4] + 1 ) < $hGeneName{$sGeneName}[2] );
    $hGeneName{$sGeneName}[3] = $aFields[5]
      if ( $hGeneName{$sGeneName}[3] < $aFields[5] );

    say STDERR "\t\t"
      . join( "\t", "Processing", $sRefSeqName, $sChrom, $sStrand, $nExonCnt )
      if $bDebug;
    #
    # Note: this is my re-work of his code
    # Overview of RefGene Processing Approach
    #
    # 1 - make arrays
    #       @this_transcript_seq            = array of bases
    #       @this_transcript_pos            = array of positions for bases
    # 2 - check for errors
    #       ensure arrays are matching lengths
    #       check whether coding transcripts start with ATG
    #       check whether coding transcripts are divisible by 3 exactly
    # 3 - save annotation for positions
    #
    #  From the mRNA:
    #     Tstart CStart        CEnd          TEnd
    #      +--------+-------------+-------------+
    #       55555555 CCCCCCCCCCCCC 3333333333333
    #
    #  5'UTR -> nPos > CStart && nPos < TStart
    #  coding -> nPos >= CStart && nPos >= CEnd
    #  3'UTR -> nPos < CEnd && nPos < TEnd
    #
    #  For reverse oriented transcripts simply reverse the transcript and swap the 3's with 5's
    #

    {
      for ( $i = 0; $i < @aExonStarts; $i++ ) {
        $nEStart = $aExonStarts[$i];
        $nEEnd   = $aExonEnds[$i] - 1;

        say STDERR "\n\t\t"
          . join( " ",
          "Processing Exon: ",
          eval( $i + 1 ),
          $nEStart, $nEEnd, $nCStart, $nCEnd )
          if $bDebug;

        #
        # build 4 arrays containing annotation information
        #   - transcript sequence
        #   - transcript BP positions
        #   - array with 5', coding base, 3' annotations
        #
        my @exon = ();
        for ( $nPos = $nEStart; $nPos <= $nEEnd; $nPos++ ) {
          push @exon,                uc( substr( $sFaSequence, $nPos, 1 ) );
          push @this_transcript_seq, uc( substr( $sFaSequence, $nPos, 1 ) );
          push @this_transcript_pos, $nPos;

          if ( $nPos < $nCEnd ) {
            if ( $nPos >= $nCStart ) {
              push @this_coding_seq, uc( substr( $sFaSequence, $nPos, 1 ) );
            }
            else {
              push @this_coding_seq, 5;
            }
          }
          else {
            push @this_coding_seq, 3;
          }
        }

        # print exon for debugging purposes
        if ($bDebug) {

          my $exon = join( "", @exon );
          if ( $sStrand eq "-" ) {
            $exon = reverse $exon;
            $exon =~ tr/ACGT/TGCA/;
          }
          say "\n\t\t>>>>" . $exon;
        }

        #
        # Annotate splice donor/acceptor bp
        #  - i.e., bp within 6 bp of exon start / stop
        #
        # From the gDNA:
        #
        #        EStart    CStart          EEnd       EStart    EEnd      EStart   CEnd      EEnd
        #        +-----------+---------------+-----------+--------+---------+--------+---------+
        #  Exons  111111111111111111111111111             22222222           333333333333333333
        #  Code               *******************************************************
        #  APR                                        ###                ###
        #  DNR                                %%%                  %%%
        #
        for ( $n = 1; $n <= SPLICE_LEN; $n++ ) {
          if ( $nEStart > $nCStart ) {
            $nSplPos = $nEStart - $n;
            if (
              !(
                ( defined $hData{$nSplPos}{ $hGeneName{$sGeneName}[1] }[1] )
                && (
                  $hData{$nSplPos}{ $hGeneName{$sGeneName}[1] }[1] == ( $sStrand eq "+" )
                  ? FWD_CODING
                  : REV_CODING
                )
              )
              )
            {

              # Assigning Gene Number
              $hData{$nSplPos}{ $hGeneName{$sGeneName}[1] }[0] = $hGeneName{$sGeneName}[0];

              # Coding/UTR & Orientation
              $hData{$nSplPos}{ $hGeneName{$sGeneName}[1] }[1] =
                ( $sStrand eq "+" ) ? FWD_SP_APR : REV_SP_APR;

              # Codon
              $hData{$nSplPos}{ $hGeneName{$sGeneName}[1] }[2] = NON_CODING;

              # Amino Acid Number
              $hData{$nSplPos}{ $hGeneName{$sGeneName}[1] }[3] = 0;

              # Error Code
              $hData{$nSplPos}{ $hGeneName{$sGeneName}[1] }[4] = $nError;

              # RefSeq Id
              $hData{$nSplPos}{ $hGeneName{$sGeneName}[1] }[5] = $sRefSeqName;

              # sha1 digest - because RefSeqID isn't always unique
              $hData{$nSplPos}{ $hGeneName{$sGeneName}[1] }[6] = $sha1_digest;

            }
          }

          if ( $nEEnd < $nCEnd ) {
            $nSplPos = $nEEnd + $n;
            if (
              !(
                ( defined $hData{$nSplPos}{ $hGeneName{$sGeneName}[1] }[1] )
                && (
                  $hData{$nSplPos}{ $hGeneName{$sGeneName}[1] }[1] == ( $sStrand eq "+" )
                  ? FWD_CODING
                  : REV_CODING
                )
              )
              )
            {

              # Assigning Gene Number
              $hData{$nSplPos}{ $hGeneName{$sGeneName}[1] }[0] = $hGeneName{$sGeneName}[0];

              # Coding/UTR & Orientation
              $hData{$nSplPos}{ $hGeneName{$sGeneName}[1] }[1] =
                ( $sStrand eq "+" ) ? FWD_SP_DNR : FWD_SP_DNR;

              # Codon
              $hData{$nSplPos}{ $hGeneName{$sGeneName}[1] }[2] = NON_CODING;

              # Amino Acid Number
              $hData{$nSplPos}{ $hGeneName{$sGeneName}[1] }[3] = 0;

              # Error Code
              $hData{$nSplPos}{ $hGeneName{$sGeneName}[1] }[4] = $nError;

              # RefSeq Id
              $hData{$nSplPos}{ $hGeneName{$sGeneName}[1] }[5] = $sRefSeqName;

              # sha1 digest - because RefSeqID isn't always unique
              $hData{$nSplPos}{ $hGeneName{$sGeneName}[1] }[6] = $sha1_digest;
            }
          }
        }
      }

      # Reverse strand adjustments
      #  - basically, get the reverse complement of the strand and swap the 5' with 3' distinction
      if ( $sStrand eq "-" ) {
        my ( @rev_this_transcript_seq, @rev_this_transcript_pos, @rev_this_coding_seq );
        my $x = 0;
        for ( my $y = $#this_transcript_seq; $y >= 0; $y-- ) {
          $rev_this_transcript_seq[$x] = $this_transcript_seq[$y];
          $rev_this_transcript_pos[$x] = $this_transcript_pos[$y];
          $rev_this_coding_seq[$x]     = $this_coding_seq[$y];
          $x++;
        }

        my $this_transcript_seq = join( "", @rev_this_transcript_seq );
        my $this_coding_seq     = join( "", @rev_this_coding_seq );
        $this_transcript_seq =~ tr/ACGT/TGCA/;
        $this_coding_seq =~ tr/ACGT53/TGCA35/;

        @this_transcript_seq = split( //, $this_transcript_seq );
        @this_transcript_pos = @rev_this_transcript_pos;
        @this_coding_seq     = split( //, $this_coding_seq );
      }
      say STDERR "Finished Generating transcript for $sStrand $sRefSeqName:\n" if $bDebug;
      say STDERR join( "\n\n",
        join( "", @this_transcript_seq ),
        join( ",", $this_transcript_pos[0], $this_transcript_pos[-1] ),
        join( "",  @this_coding_seq ) )
        if $bDebug;

      # check for fatal errors
      die unless ( scalar @this_transcript_pos == scalar @this_transcript_seq );

      # check for soft errors in transcript
      #  1 : out of frame
      #  2 : doesn't start with ATG
      #  3 : if both are the case
      if ( abs( $nCStart - $nCEnd ) > 0 ) {
        my $seq = join( "", @this_coding_seq );
        $seq =~ s/3|5//g; # remove the 5' or 3' distinction

        if ( !( ( length $seq ) % 3 == 0 ) ) {
          $nError += 1;
          warn "$sStrand $sRefSeqName: transcript not divisible by 3\n\n>$seq\n\n";
        }
        if ( $seq !~ m/^ATG/ ) {
          $nError += 2;
          warn "$sStrand $sRefSeqName ($sha1_digest): does not start with ATG...\n"
            . join( "\t", @aFields ) . "\n"
            . "\n\n>$seq\n\n";
        }
      }
      else {
        $nError = 0;
      }

      # process transcripts
      my $coding_bp = 0;
      foreach my $j ( 0 .. $#this_transcript_pos ) {
        $nPos = $this_transcript_pos[$j] + 1;
        my $site_code = uc $this_coding_seq[$j];

        # Assigning Gene Number
        $hData{$nPos}{ $hGeneName{$sGeneName}[1] }[0] = $hGeneName{$sGeneName}[0];

        # Coding/UTR & Orientation && Codon && Amino Acid Number
        if ( $site_code eq "5" ) {
          $hData{$nPos}{ $hGeneName{$sGeneName}[1] }[1] =
            ( $sStrand eq "+" ) ? FWD_5P_UTR : REV_5P_UTR;
          $hData{$nPos}{ $hGeneName{$sGeneName}[1] }[2] = NON_CODING;
          $hData{$nPos}{ $hGeneName{$sGeneName}[1] }[3] = 0;
        }
        elsif ( $site_code eq "3" ) {
          $hData{$nPos}{ $hGeneName{$sGeneName}[1] }[1] =
            ( $sStrand eq "+" ) ? FWD_3P_UTR : REV_3P_UTR;
          $hData{$nPos}{ $hGeneName{$sGeneName}[1] }[2] = NON_CODING;
          $hData{$nPos}{ $hGeneName{$sGeneName}[1] }[3] = 0;
        }

        elsif ( $site_code =~ m/[ACGT]/ ) {
          $hData{$nPos}{ $hGeneName{$sGeneName}[1] }[1] =
            ( $sStrand eq "+" ) ? FWD_CODING : REV_CODING;

          # note that if there is an error in the codon then the nError will already be set from above nError = 1 -> improper length / out of frame
          my ( $codon_num, $codon_pos, $codon_start, $codon_stop, $codon_seq, $codon_code,
            $codon_seq_pos_code );
          $codon_num   = int( ( $coding_bp / 3 ) );
          $codon_pos   = eval( $coding_bp % 3 );
          $codon_start = $j - $codon_pos;
          $codon_stop  = $codon_start + 2;
          if ( defined $this_transcript_seq[$codon_start]
            && defined $this_transcript_seq[$codon_stop] )
          {
            $codon_seq = uc( join( "", @this_transcript_seq[ $codon_start .. $codon_stop ] ) );
            $codon_code = Codon_Code( \%hCmdLineOption, $codon_seq );

            # Combined codon code and the code for the BP position within the codon
            $codon_seq_pos_code = $codon_code + $codon_pos * 64;

            $hData{$nPos}{ $hGeneName{$sGeneName}[1] }[2] = $codon_seq_pos_code;
            $hData{$nPos}{ $hGeneName{$sGeneName}[1] }[3] = $codon_num + 1;
          }
          else {
            $hData{$nPos}{ $hGeneName{$sGeneName}[1] }[2] = NON_CODING;
            $hData{$nPos}{ $hGeneName{$sGeneName}[1] }[3] = 0;
          }
          $coding_bp++;

          print STDERR join( "\t", $j, $nPos, $site_code, $nError, $sRefSeqName ) if $bDebug;
          print STDERR "\t"
            . join( "\t",
            $codon_num, $codon_pos, $codon_start, $codon_stop, $codon_seq, $codon_code )
            if $bDebug;
          print STDERR "\n" if $bDebug;
        }
        else {
          die "Something bad happened with the annotation...\n";
        }

        # Error Code
        $hData{$nPos}{ $hGeneName{$sGeneName}[1] }[4] = $nError;

        # RefSeq Id
        $hData{$nPos}{ $hGeneName{$sGeneName}[1] }[5] = $sRefSeqName;

        # sha1 digest - because RefSeqID isn't always unique
        $hData{$nPos}{ $hGeneName{$sGeneName}[1] }[6] = $sha1_digest;

        #print STDERR "\n";
      }
    }
    $nX++;

    # for testing purposes
    #die if ($sRefSeqName eq "NM_001171136");
    #die if ($sRefSeqName eq "NM_032971");
    #die if ($sRefSeqName eq "NM_001291462");
  }
  close($fh);
}

sub BuildBinary {

  $fpGID = $sGeneFile->filehandle(">");
  binmode $fpGID;

  my $xGene;
  for ( $i = 0; $i < @aGenes; $i++ ) {
    $sGeneName = $aGenes[$i];
    if ( !defined $hGeneCoords{ $hGeneName{$sGeneName}[2] } ) {
      $hGeneCoords{ $hGeneName{$sGeneName}[2] } = $hGeneName{$sGeneName}[3];
    }

    $hGeneCoords{ $hGeneName{$sGeneName}[2] } = (
      ( $hGeneCoords{ $hGeneName{$sGeneName}[2] } < $hGeneName{$sGeneName}[3] )
      ? $hGeneName{$sGeneName}[3]
      : $hGeneCoords{ $hGeneName{$sGeneName}[2] }
    );

    $xGene = pack( "Z120", $sGeneName );
    print $fpGID "$xGene";
  }

  close($fpGID);

  $fpOUT = $sOutFile->filehandle(">");
  $fpSNP = $sSnpFile->filehandle(">");
  $fpIDX = $sIdxFile->filehandle(">");

  binmode $fpOUT;
  binmode $fpSNP;
  binmode $fpIDX;

  my $xRecord;
  my ( @aRsNum,      @aFreq,         @aOrientation );
  my ( $nOutFilePos, $nOutRecordLen, $nSnpFilePos, $nSnpRecordLen );
  my ( $nGeneStart,  $nGeneEnd,      $nNonCoding );
  my $sTypeDef;
  my @aGeneStarts = sort { $a <=> $b } keys %hGeneCoords;
  my $nI = 0;

  $nOutFilePos = $nSnpFilePos = 0;

  $nNonCoding = INTERGENIC;

  if (@aGeneStarts) {
    $nGeneStart = $aGeneStarts[$nI];
    $nGeneEnd   = $hGeneCoords{$nGeneStart};
  }
  else {
    $nGeneStart = $nGeneEnd = $nFaSeqLen + 100;
  }

  for ( $nPos = 1; $nPos <= $nFaSeqLen; $nPos++ ) {

    if ( $nGeneStart <= $nPos ) {
      $nNonCoding = INTRONIC;
    }

    $nOutRecordLen = $nSnpRecordLen = 0;

    if ( defined $hData{$nPos} ) {
      foreach my $nGene ( sort { $a <=> $b } keys %{ $hData{$nPos} } ) {
        $sTypeDef = "I C C S C Z25 h40";

        $nOutRecordLen += length( pack( $sTypeDef, () ) );

        $xRecord = pack( $sTypeDef,
          $hData{$nPos}{$nGene}[0], $hData{$nPos}{$nGene}[1], $hData{$nPos}{$nGene}[2],
          $hData{$nPos}{$nGene}[3], $hData{$nPos}{$nGene}[4], $hData{$nPos}{$nGene}[5],
          $hData{$nPos}{$nGene}[6] );

        print $fpOUT "$xRecord";
      }

      $nNonCoding = CODING_POS;
    }

    if ( ( defined $hCmdLineOption{'dbsnpdir'} )
      && ( defined $hCmdLineOption{'snpbld'} ) )
    {
      if ( defined $hdbSnp{$nPos} ) {
        $sTypeDef = "L C C";

        @aRsNum       = split( /,/, $hdbSnp{$nPos}[0] );
        @aFreq        = split( /,/, $hdbSnp{$nPos}[1] );
        @aOrientation = split( /,/, $hdbSnp{$nPos}[2] );

        for ( $i = 0; $i < @aRsNum; $i++ ) {
          $nSnpRecordLen += length( pack( $sTypeDef, () ) );

          $xRecord = pack( $sTypeDef, $aRsNum[$i], $aFreq[$i], $aOrientation[$i] );

          print $fpSNP "$xRecord";
        }
      }
    }

    if ( $nOutRecordLen >= 65_535 ) {
      die "$nOutRecordLen is too large for a short ...\n";
    }

    $sTypeDef = "L S L S C";

    $xRecord = pack( $sTypeDef,
      $nOutFilePos, $nOutRecordLen, $nSnpFilePos, $nSnpRecordLen, $nNonCoding );

    print $fpIDX "$xRecord";

    $nOutFilePos += $nOutRecordLen;
    $nSnpFilePos += $nSnpRecordLen;

    if ( $nPos == $nGeneEnd ) {
      $nNonCoding = INTERGENIC;

      while ( $nGeneEnd <= $nPos ) {
        $nI++;
        if ( $nI < @aGeneStarts ) {
          $nGeneStart = $aGeneStarts[$nI];
          $nGeneEnd   = $hGeneCoords{$nGeneStart};
        }
        else {
          $nGeneStart = $nGeneEnd = $nFaSeqLen + 100;
        }
      }
    }
  }

  close($fpIDX);
  close($fpSNP);
  close($fpOUT);
}

# Generate_Codon_Pos_Code()
#
# Purpose
#	generates code for the codon and the position within the codon of
#	the base coordinate.
#
# Required Parameters
#   phCmdLineOption	= pointer to hash containing command line options
#	sSequence		= chromosome sequence
#   nStart			= coding start position
#   nEnd			= coding stop position
#	sStrand			= exon orientation
#	nBasePos		= base position
#
# Optional Parameters
#   none
#
# Returns
#   nCodonCode		= code for codon and base position within codon
#
# Side Effects
#   none
#
# Assumptions
#
# Notes
#

sub Generate_Codon_Pos_Code {
  my $phCmdLineOption = shift;
  my $sSequence       = uc(shift);
  my $nStart          = shift;
  my $nEnd            = shift;
  my $sStrand         = shift;
  my $nBasePos        = shift;

  my $sSubName = ( caller(0) )[3];

  if (
    !(
         ( defined $sSequence )
      && ( defined $nStart )
      && ( defined $nEnd )
      && ( defined $sStrand )
      && ( defined $nBasePos )
      && ( defined $phCmdLineOption )
    )
    )
  {
    croak "$sSubName - required parameters missing\n";
  }

  # Local variables
  my $bDebug   = ( defined $phCmdLineOption->{'debug'} )   ? TRUE : FALSE;
  my $bVerbose = ( defined $phCmdLineOption->{'verbose'} ) ? TRUE : FALSE;

  my ( $nCodonPos, $nCodonStart, $sCodon, $nCodonCode );

  # Start

  if ( $sStrand =~ /^\+$/ ) {
    $nCodonPos = ( $nBasePos - $nStart + 1 ) % 3;
    if ( $nCodonPos == 0 ) { $nCodonPos = 3; }

    $nCodonStart = ( ( int( ( ( $nBasePos - $nStart + 1 ) / 3 ) + 0.70 ) ) - 1 ) * 3;

    $sCodon = substr( $sSequence, $nCodonStart, 3 );
  }
  else {
    $nCodonPos = ( $nStart - $nBasePos + 1 ) % 3;
    if ( $nCodonPos == 0 ) { $nCodonPos = 3; }

    $nCodonStart = -1 * ( ( int( ( ( $nStart - $nBasePos + 1 ) / 3 ) + 0.70 ) ) * 3 );

    $sCodon = substr( $sSequence, $nCodonStart, 3 );
    $sCodon = reverse($sCodon);
    $sCodon =~ tr/ACGT/TGCA/;
  }

  $nCodonCode = Codon_Code( $phCmdLineOption, $sCodon );
  $nCodonCode += ( ( $nCodonPos - 1 ) * 64 );

  return $nCodonCode;
}

# Codon_Code()
#
# Purpose
#   generates a number between 1 and 64 based on codon.
#
# Required Parameters
#   phCmdLineOption = pointer to hash containing command line options
#   sCodon			= codon sequence
#
# Optional Parameters
#   none
#
# Returns
#   nCodonCode		= code for codon and base position within codon
#
# Side Effects
#   none
#
# Assumptions
#
# Notes
#
sub Codon_Code {
  my $phCmdLineOption = shift;
  my $sCodon          = shift;

  my $sSubName = ( caller(0) )[3];

  if ( !( ( defined $sCodon ) && ( defined $phCmdLineOption ) ) ) {
    croak "$sSubName - required parameters missing\n";
  }

  # Local variables
  my $bDebug   = ( defined $phCmdLineOption->{'debug'} )   ? TRUE : FALSE;
  my $bVerbose = ( defined $phCmdLineOption->{'verbose'} ) ? TRUE : FALSE;

  my ( @aCodon, $nCodonCode );

  # Start

  if ( $sCodon =~ /N/ ) {
    return 0;
  }
  else {
    $sCodon =~ tr/ACGT/0123/;
    @aCodon = split( //, $sCodon );
    $nCodonCode = 0;

    for ( my $c = 0; $c < @aCodon; $c++ ) {
      $nCodonCode += ( $aCodon[$c] * ( 4**( 2 - $c ) ) );
    }
    $nCodonCode += 1;

    return $nCodonCode;
  }
}

# Parse_dbSnp()
#
# Purpose
#   parses the DbSnp flat file for SNP information
#
# Required Parameters
#   phCmdLineOption = pointer to hash containing command line options
#   phSnp			= pointer to hash to store SNP information
#
# Optional Parameters
#   none
#
# Returns
#   none
#
# Side Effects
#   modifies the hash containing SNP information
#
# Assumptions
#
# Notes
#
sub Parse_dbSnp {
  my $phCmdLineOption = shift;
  my $phSnp           = shift;

  my $sSubName = ( caller(0) )[3];
  if ( !( ( defined $phCmdLineOption ) && ( defined $phSnp ) ) ) {
    croak "$sSubName - required parameters missing\n";
  }

  # Local variables
  my $bDebug   = ( defined $phCmdLineOption->{'debug'} )   ? TRUE : FALSE;
  my $bVerbose = ( defined $phCmdLineOption->{'verbose'} ) ? TRUE : FALSE;

  my $sGenome    = $phCmdLineOption->{'genome'};
  my $sChr       = "chr" . $phCmdLineOption->{'chr'};
  my $sSnpBuild  = $phCmdLineOption->{'snpbld'};
  my $sDbSnpFile = path( $phCmdLineOption->{'dbsnpdir'} )
    ->child("$sGenome.$sChr.dbsnp$sSnpBuild.txt");
  my $fh = $sDbSnpFile->filehandle("<");

  my $nX       = 0;
  my $nSnpId   = 0;
  my $nHetRate = 0;
  my $nO       = 0;
  while (<$fh>) {
    $_ =~ s/\s+$//g; # remove all spaces
    next if ( $_ =~ /^#/ ); # skip comments

    # this is needlessly complex
    my ( $sChrom, $nSnpStart, $nSnpEnd, $sSnpId, $sSnpStrand, $sSnpAlleleFreqs ) =
      ( split( /\t/, $_ ) )[ 1, 2, 3, 4, 6, 24 ];

    my $maf_allele_freq = 0;

    if ( $sSnpId =~ /^rs(\d+)/ ) {
      $nSnpId = $1;

      foreach my $nP ( ( $nSnpStart + 1 ) .. ($nSnpEnd) ) {
        if ($sSnpAlleleFreqs) {
          my @subfields = split( /,/, $sSnpAlleleFreqs );
          my @aSnpAlleleFreqs = sort { $b <=> $a } @subfields;
          $maf_allele_freq = sprintf( "%0.6f", eval( 1 - $aSnpAlleleFreqs[0] ) );
        }

        my $maf_code = int( $maf_allele_freq * 199 ) + 1;
        if ( $sSnpStrand =~ /^\+$/ ) {
          $nO = 1;
        }
        elsif ( $sSnpStrand =~ /^\-$/ ) {
          $nO = 2;
        }
        else {
          $nO = 0;
        }

        if ( $sChrom =~ /$sChr/ ) {
          $phSnp->{$nP}[0] .= "$nSnpId,";
          $phSnp->{$nP}[1] .= "$maf_code,";
          $phSnp->{$nP}[2] .= "$nO,";

          #say STDERR join("\t", $sChr, $nP, @{ $phSnp->{$nP}});
        }
      }
    }
    $nX++;
  }
  close($fh);
}

sub Print_Sequence {
  my $sSeq = uc(shift);
  my $nBrk = shift;

  my $sSubSeq = "";
  for ( my $s = 0; $s < length($sSeq); $s += $nBrk ) {
    if ( ( $s + $nBrk ) > length($sSeq) ) {
      $sSubSeq = substr( $sSeq, $s, ( length($sSeq) - $s ) );
    }
    else {
      $sSubSeq = substr( $sSeq, $s, $nBrk );
    }
    print STDERR "\t\t$sSubSeq\n";
  }
}

# CleanChr takes a string, assumed to be formatted like chrX, Chr1, 15, or 24
# and returns 'chr1' .. 'chr22', 'chrX', etc.
# Using plink as a guide for num -> letter chromosomes:
# http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml
# The autosomes should be coded 1 through 22. The following other codes can be
# used to specify other chromosome types:
#    X    X chromosome                    -> 23
#    Y    Y chromosome                    -> 24
#    XY   Pseudo-autosomal region of X    -> 25
#    MT   Mitochondrial                   -> 26
sub CleanChr {
  my $chr = shift;

  my @chrs = ( 1 .. 22, 'X', 'Y', 'M' );
  my %chrOk = map { $_ => 1 } (@chrs);
  my %chrNumToLetter = ( 23 => 'X', 24 => 'Y', 26 => 'M' );

  $chr = uc($chr);
  $chr =~ s/\ACHR//xm;
  if ( exists $chrOk{$chr} ) {
    return "chr" . $chr;
  }
  elsif ( exists $chrNumToLetter{$chr} ) {
    return "chr" . $chrNumToLetter{$chr};
  }
  else {
    my $msg = "Error: unrecognized chromosome: $chr";
    say $msg;
    exit(1);
  }
}

##############################################################################
### POD Documentation
##############################################################################

__END__

=head1 NAME

create_db.pl

=head1 SYNOPSIS

    create_db.pl
      --g <genome_build>
      --c <chromosome>
      --f <dir with fasta files>
      --genedir <dir with gene files>
      --geneset <refGene or knownGene>
      --o <output_directory>
      [--d <dir with snp files>]
      [--n <SNP_build>]
      [--v]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

    --g <genome_build>                = genome whose refgene file exists.
    --c <chromosome>                  = chromosome number.
    --f <dir with fasta files>        = path to the directory containing genome directory of fasta files.
    --genedir <dir with gene files>   = path to the directory containing genome directory of refgene files.
    --geneset <refGene or knownGene>  = type of gene data
    --d <dir with dbsnp files>        = path to the directory containing genome directory of dbsnp files.
    --n <SNP_build>                   = dbSNP build.
    --o <output_dir>                  = output directory.
    --v                               = generate runtime messages. Optional

=head1 DESCRIPTION

This program creates a binary representation of the genome that incorporates gene
and snp data. The goal was to modify Amol's original code as little as possible.

=head1 DEPENDENCIES

Carp
Getopt::Long
Pod::Usage
Digest::SHA
Path::Tiny

=head1 BUGS AND LIMITATIONS

All software has bugs. If you find one please email Thomas Wingo
(<thomas.wingo@emory.edu>).

=head1 AUTHOR

Thomas S. Wingo based on code written by Amol Shetty.

=head1 LICENSE AND COPYRIGHT

Copyright (c) 2016 Thomas S. Wingo (<thomas.wingo@emory.edu>). All rights
reserved.

This program is free software; you can distribute it and/or modify it under the
same terms as GNU GPL 3.0.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or FITNESS
FOR A PARTICULAR PURPOSE. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=cut
