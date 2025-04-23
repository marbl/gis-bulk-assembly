
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
use Cwd 'abs_path';

@ISA    = qw(Exporter);
@EXPORT = qw(getParameters loadSamples numFiles dataAvailable $root $rasm $rsoft $refn $refc $odb %samples);

use strict;
use warnings "all";
no  warnings "uninitialized";

#  These are also set in analyze.sh.

our $root  = abs_path('.');    #"/data/Phillippy2/projects/hprc-assemblies";
our $rasm  = "/invalid/path";  # now set in main; "/data/Phillippy2/projects/hprc-assemblies/assemblies-v3";
our $rsoft = "/invalid/path";  # now set tin main;
our $refn  = "/invalid/path";  # the reference non-HPC
our $refc  = "/invalid/path";  # the reference HPC
our $odb   = "/invalid/path";  # the ODB database to use (e.g. primates_odb10)
our %samples;

$ENV{'REF_CACHE'} = "$root/hprc-cache/samtools";
$ENV{'HPRC_ROOT'} = "$root";

#
#  Extract the values from within brackets: "[ A,B,C,D ]"
#  Returns reference to an array.
#  If given a non-bracketed string, returns the original string.
#
sub extractItems ($) {
  my $orig = shift @_;

  return undef   if ($orig eq "na");
  return undef   if ($orig eq "NA");
  return undef   if ($orig eq "");
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

      $h =~ s/parameters/params/;
      $h =~ s/options/params/;

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

#
#  Return optional parameters if the sample has them defined
#  Otherwise return empty string
#
sub getParameters ($) {
  my $samp  = shift @_;

  my $opts  = $samples{$samp}{"params"};
  return (defined($opts)) ? join(" ", @$opts) : ""
}

#
#  Return the number of files for a given sample and type.
#  This is NOT the number of files that are present locally;
#  use getFileMap() to get that list.
#
sub numFiles ($$) {
  my $samp  = shift @_;
  my $type  = shift @_;

  $type = "hifi"   if ($type eq "hifi-cutadapt");

  my $fles  = $samples{$samp}{$type};
  return (defined($fles)) ? scalar(@$fles) : 0
}




sub dataAvailable ($@) {
  my $samp = shift @_;
  my $flav = shift @_;

  my $tu  = ((numFiles($samp, "mat-ilmn") == 0) || (numFiles($samp, "pat-ilmn") == 0));
  my $hu  =  (numFiles($samp, "hic")      == 0);

  my $uT = "";
  my $uH = "";
  my $uC = "";

  if    ($tu && $hu) {
    $uT = "no trio data";
    $uH = "no hi-c data";
    $uC = "no trio or hi-c data";
  }
  elsif ($tu) {
    $uT = "no trio data";
    $uC = "no trio data";
  }
  elsif ($hu) {
    $uH = "no hi-c data";
    $uC = "no hi-c data";
  }

  if (!defined($flav))           { return($uT, $uH, $uC); }
  elsif ($flav eq "verkko-trio") { return $uT; }
  elsif ($flav eq "verkko-hi-c") { return $uH; }
  elsif ($flav eq "verkko-thic") { return $uC; }
  else                           { return "";  }
}

1;
