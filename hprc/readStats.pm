
###############################################################################
#
#  This file is part of a pipeline to (help) run Verkko assemblies on
#  numerous HPRC samples.
#
#  This is a 'United States Government Work', and is released in the public
#  domain.
#
##

package hprc::readStats;
require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(summarizeReads);

use strict;
use warnings "all";
no  warnings "uninitialized";

use hprc::samples;
use hprc::aws;

sub summarizeReads ($$$) {
  my $s      = shift @_;
  my $types  = shift @_;
  my $opts   = shift @_;

}

1;
