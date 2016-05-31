use 5.10.0;
use strict;
use warnings;

package GenPro;

# ABSTRACT: Make personal protein databases using genomic variant information.

our $VERSION = '0.01';

# Dependencies
use strict;
use warnings;
use base qw/ Exporter /;
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

#use Type::Params qw/ compile /;
#use Types::Standard qw/ Str ArrayRef HashRef RegexpRef Int Num /;

our @EXPORT    = qw/ /;
our @EXPORT_OK = qw/ /;
