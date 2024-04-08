
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
#  If given a non-bracketed string, returns the original string.
#
sub extractItems ($) {
  my $orig = shift @_;

  return undef   if ($orig eq "na");
  return $orig   if ($orig !~ m/^\s*\[/);

  $orig =~ s/^\s*\[\s*//;
  $orig =~ s/\s*\]\s*$//;

  my $items = [];
  foreach my $item (split ',\s*', $orig) {
    $item =~ s/^\s*['`"]//;
    $item =~ s/['`"]\s*$//;

    push @$items, $item;
  }

  return $items;
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
  while (!eof(F)) {
    $_ = <F>;
    s/^\s+//;
    s/\s+$//;

    next if (($_ eq "") || ($_ =~ m/^#/));

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

    last;
  }

  while (!eof(F)) {
    $_ = <F>;
    s/^\s+//;
    s/\s+$//;

    next if (($_ eq "") || ($_ =~ m/^#/));

    my @vals = split '\t', $_;
    my $s = {};

    my $col = 0;
    my $id;

    if (scalar(@vals) != scalar(@cols)) {
      print "ERROR: Expecting cols=", scalar(@cols), " values, but got ", scalar(@vals), ".\n";
      print "ERROR: $_\n";
    }

    while (scalar(@vals) > 0) {
      $$s{$cols[$col++]} = extractItems(shift @vals);
    }

    $id = $$s{'sample_id'};

    die if (scalar(@vals) > 0);

    $samples{$id} = $s;
  }
}

1;
