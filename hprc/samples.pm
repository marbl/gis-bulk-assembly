
###############################################################################
#
#  This file is part of a pipeline to (help) run Verkko assemblies on
#  numerous HPRC samples.
#
#  This is a 'United States Government Work', and is released in the public
#  domain.
#
##

package hprc::samples;
require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(loadSamples $root %samples);

use strict;
use warnings "all";
no  warnings "uninitialized";

our $root = -e "/data/walenzbp/hprc" ? "/data/walenzbp/hprc" : "/work/hprc";
our %samples;

#
#  Extract the values from within brackets: "[ A,B,C,D ]"
#  Returns reference to an array.
#
sub extractItems ($) {
  my $orig = shift @_;

  $orig =~ s/^\s*\[\s*//;
  $orig =~ s/\s*\]\s*$//;

  my $items = [];
  foreach my $item (split ',\s*', $orig) {
    $item =~ s/^\s*['`"]//;
    $item =~ s/['`"]\s*$//;

    push @$items, $item;
  }

  return($items);
}

#
#  Read a .tsv file containing sample information.
#
sub loadSamples ($) {
  my $tsv = shift @_;
  my @header;
  my @cols;

  open(F, "< $tsv") or die "Failed to open '$tsv': $!\n";

  #  Read the header line and massage into more convenient labels.
  {
    $_ = <F>;
    s/^\s+//;
    s/\s+$//;

    foreach my $h (split '\t', $_) {
      my $o = $h;

      $h =~ tr/A-Z/a-z/;
      $h =~ s/\s+/_/g;

      $h =~ s/paternal/pat/;
      $h =~ s/maternal/mat/;

      $h =~ s/population/pop/;
      $h =~ s/hifi_product_site/prod_site/;
      $h =~ s/production_year/prod_year/;

      $h =~ s/ontcov.over100kb/ont_cov_100k/;

      $h =~ s/nanopore/ont/;

      $h =~ s/child_ilmn/ilmn/;
      $h =~ s/mat_ilmn/mat-ilmn/;
      $h =~ s/pat_ilmn/pat-ilmn/;

      #printf "%20s -> %s\n", $o, $h;

      push @cols, $h;
    }
  }

  while (<F>) {
    my @vals = split '\t', $_;
    my $s = {};

    my $col = 0;
    my $id;

    $$s{$cols[$col++]} = $id        = shift @vals;
    $$s{$cols[$col++]}              = shift @vals;
    $$s{$cols[$col++]}              = shift @vals;
    $$s{$cols[$col++]}              = shift @vals;
    $$s{$cols[$col++]}              = shift @vals;
    $$s{$cols[$col++]}              = shift @vals;
    $$s{$cols[$col++]}              = shift @vals;
    $$s{$cols[$col++]}              = shift @vals;
    $$s{$cols[$col++]}              = shift @vals;
    $$s{$cols[$col++]}              = shift @vals;
    $$s{$cols[$col++]}              = shift @vals;
    $$s{$cols[$col++]} = extractItems(shift @vals);
    $$s{$cols[$col++]} = extractItems(shift @vals);
    $$s{$cols[$col++]} = extractItems(shift @vals);
    $$s{$cols[$col++]} = extractItems(shift @vals);
    $$s{$cols[$col++]} = extractItems(shift @vals);
    $$s{$cols[$col++]} = extractItems(shift @vals);

    die if (scalar(@vals) > 0);

    $samples{$id} = $s;
  }
}

1;
