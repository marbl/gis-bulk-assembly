#!/usr/bin/env perl

###############################################################################
#
#  This file is part of a pipeline to (help) run Verkko assemblies on
#  numerous HPRC samples.
#
#  This is a 'United States Government Work', and is released in the public
#  domain.
#
##

use strict;
use warnings "all";
no  warnings "uninitialized";

use FindBin;
use lib "$FindBin::RealBin";

use hprc::samples;
use hprc::aws;
use hprc::list;
use hprc::readStats;
use hprc::readFilter;
use hprc::status;
use hprc::hapmers;
use hprc::assemble;


my @errs;
my $mode = (scalar(@ARGV) > 0) ? shift @ARGV : "help";
my %opts;

my %sampleList;
my %readTypes;

#
#  Load the sample list.
#

#oadSamples("batch-1-test.tsv");
loadSamples("batch-2-with-hic.tsv");

#
#  Parse args.
#

while (scalar(@ARGV) > 0) {
  my $arg = $ARGV[0];  shift @ARGV;

  if ($arg eq "--sample")  {
    while ((scalar(@ARGV) > 0) &&      #  While more words
           ($ARGV[0] !~ m/^-/)) {      #  and not an option word,
      $sampleList{ shift @ARGV } = 1;  #  add sample to the list
    }
  }

  elsif ((($mode eq "fetch") ||
          ($mode eq "list")) && ($arg eq "--type")) {
    while ((scalar(@ARGV) > 0) &&     #  While more words
           ($ARGV[0] !~ m/^-/)) {     #  and not an option word,
      $readTypes{ shift @ARGV } = 1;  #  add data-type to the list
    }
  }

  elsif (($mode eq "list") && ($arg eq "--files")) {
    $opts{"files"} = 1;
  }

  elsif (($mode eq "assemble") && ($arg eq "--canu-trio"))  { $opts{"canu-trio"}   = 1; }   #  The default is to run both
  elsif (($mode eq "assemble") && ($arg eq "--canu-hifi"))  { $opts{"canu-hifi"}   = 1; }   #  verkko-trio then verkko-hi-c.
  elsif (($mode eq "assemble") && ($arg eq "--trio"))       { $opts{"verkko-trio"} = 1; }   #
  elsif (($mode eq "assemble") && ($arg eq "--hi-c"))       { $opts{"verkko-hi-c"} = 1; }   #

  elsif ((($mode eq "hapmers") || ($mode eq "assemble")) && ($arg eq "--submit")) {
    $opts{"submit"} = 1;
  }

  else {
    push @errs, "unknown option '$arg'";
  }
}



#
#  Expand sample 'all' to include all known samples and check for invalid ones.
#
if ((scalar(keys %sampleList) == 0) ||   #  If no samples supplied, or
    (exists($sampleList{'all'}))) {      #  'all' explicitly requested,
  foreach my $s (keys %samples) {        #  add all known samples to our list.
    $sampleList{$s} = 1;
  }

  delete $sampleList{'all'};
}
{
  my @missing;

  foreach my $s (sort keys %sampleList) {
    push @missing, $s   if (!exists($samples{$s}));
  }

  if (scalar(@missing) > 0) {
    my $ss = (scalar(@missing) == 1) ? "" : "s";
    my $ms = "'" . (shift @missing) . "'";
    my $ls = pop @missing;
    foreach my $s (@missing) {
      $ms .= ", '$s'";
    }
    $ms .= " and '$ls'"  if (defined($ls));

    push @errs, "Sample$ss $ms not known.\n";
  }
}



#
#  Expand fetch options and check for invalid ones.
#
if (scalar(keys %readTypes) == 0) {
  $readTypes{'all'} = 1;
}
if (exists($readTypes{'all'})) {
  $readTypes{$_} = 1   foreach qw(hifi ont hic ilmn mat-ilmn pat-ilmn);
  delete $readTypes{'all'};
}
foreach my $f (keys %readTypes) {   #  Check for invalid fetch options.
  if (($f ne 'hifi') && ($f ne 'ont') && ($f ne 'hic') && ($f ne 'ilmn') &&
      ($f ne 'mat-ilmn') && ($f ne 'pat-ilmn')) {
    delete $readTypes{$_}  foreach qw(hifi ont hic ilmn mat-ilmn pat-ilmn);
    push @errs, "Invalid '$mode' types '" . join("', '", sort keys %readTypes) . "'.";
    push @errs, "  Must be one of:";
    push @errs, "    'all'                         (everything)";
    push @errs, "    'hifi', 'ont', 'hic', 'ilmn', (child)";
    push @errs, "    'mat-ilmn', 'pat-ilmn'        (parental)";
    last;
  }
}

#
#  Check for invalid modes.
#

if (($mode ne "help") &&
    ($mode ne "list") &&
    ($mode ne "status") &&
    ($mode ne "fetch") &&
    ($mode ne "read-stats") &&
    ($mode ne "read-filter") &&
    ($mode ne "hapmers") &&
    ($mode ne "assemble") &&
    ($mode ne "postprocess")) {
  push @errs, "Invalid mode '$mode'.\n";
}

if (($mode eq "help") || (scalar(@errs) > 0)) {
  print "usage: $0 mode [options]";
  print "  MODES:\n";
  print "    list [--type all|hifi|etc] [--files]\n";
  print "    status\n";
  print "    fetch [--type all|hifi|etc]\n";
  print "    stats\n";    #  --reads or --assembly
  print "    hapmers\n";
  print "    assemble [--trio, --hi-c --canu-trio --canu-hifi] [--submit]\n";
  print "\n";
  print "  OPTIONS:\n";
  print "  --sample [sample ...]   - restrict operation to the specified samples\n";
  print "\n";
  print "ERROR: $_\n"  foreach (@errs);
  exit(1);
}


if ($mode eq "list") {
  foreach my $s (sort keys %sampleList) {
    displaySample($s, \%readTypes, \%opts);
  }
}

elsif ($mode eq "status") {
  foreach my $s (sort keys %sampleList) {
    checkStatus($s, \%opts);
  }
}

elsif ($mode eq "fetch") {
  foreach my $s (sort keys %sampleList) {
    foreach my $f (sort keys %readTypes) {
      fetchData($s, $f, \%opts);
    }
  }
}

elsif ($mode eq "read-stats") {
}

elsif ($mode eq "status") {
}

elsif ($mode eq "hapmers") {
  foreach my $s (sort keys %sampleList) {
    computeHapmers($s, ,\%opts);
  }
}

elsif ($mode eq "assemble") {
  foreach my $s (sort keys %sampleList) {
    startAssembly($s, \%opts);
  }
}

else {
}
