use 5.10.0;
use strict;
use warnings;

package GenPro;

# ABSTRACT: Make personal protein databases using genomic variant information.

our $VERSION = '0.01';

# Dependencies
use strict;
use warnings;
use Exporter 5.57 (qw/import/);
use Carp;
use Cpanel::JSON::XS;
use DB_File;
use Data::Dump;
use Getopt::Long;
use Fcntl qw/ :DEFAULT :seek /;
use Path::Tiny;
use Digest::SHA qw/ sha1_hex /;
use Getopt::Long;
use IO::Uncompress::Gunzip qw/ $GunzipError /;
use Pod::Usage;
use Path::Tiny;
use Scalar::Util qw/reftype/;

#use Type::Params qw/ compile /;
#use Types::Standard qw/ Str ArrayRef HashRef RegexpRef Int Num /;

our @EXPORT    = qw/ /;
our @EXPORT_OK = qw/ /;

# Package variables
my ( %hIUPAC, %codon_2_aa );

#<<< No perltidy
%codon_2_aa = (
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
