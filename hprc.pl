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

use v5.34;
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
use hprc::genomescope;
use hprc::hapmers;
use hprc::yakmers;
use hprc::assemble;
use hprc::analyze;
use hprc::ribotin;


my @errs;
my $mode = (scalar(@ARGV) > 0) ? shift @ARGV : "help";
my %opts;

my %sampleList;
my %readTypes;

#
#  Parse args.
#   - the first arg ('mode') is set above.
#   - the rest of the args are conditional on the mode selected,
#     though some are options are available for every mode.
#   - choice of samples and output directory is conditional
#     on options
#       --v3
#       --v4
#  All software will be unversioned and pulled from the folder named $root/software
#      but verkko and graphaligner will be pulled from the one specified by the version below

while (scalar(@ARGV) > 0) {
  my $arg = $ARGV[0];  shift @ARGV;

  #  Select a set of samples to use.
  #
  if ($arg =~ /^--v(\d+)$/) {
    my $v = $1;

    if ($v == 3) {
      loadSamples("$root/hprc-cache/b2.tsv");
      loadSamples("$root/hprc-cache/b3.tsv");
      $rasm  = "$root/assemblies-v3";
      $rsoft = "software-v3";
    }
    elsif ($v == 4) {
      loadSamples("$root/hprc-cache/b2.tsv");
      loadSamples("$root/hprc-cache/b3.tsv");
      loadSamples("$root/hprc-cache/b8.tsv");
      loadSamples("$root/hprc-cache/b1.tsv");
      loadSamples("$root/hprc-cache/b5.tsv");
      loadSamples("$root/hprc-cache/b6.tsv");
      loadSamples("$root/hprc-cache/b7.tsv");
      loadSamples("$root/hprc-cache/b9.tsv");
      loadSamples("$root/hprc-cache/b10.tsv");
      loadSamples("$root/hprc-cache/b11.tsv");
     # loadSamples("$root/hprc-cache/yr4_30hr.tsv");
      $rasm  = "$root/assemblies-v4";
      $rsoft = "software-v4";
    }
    else {
    }
  }

  #  --samples is available for every mode.
  #
  elsif ($arg =~ m/^--samples?$/)  {
    while ((scalar(@ARGV) > 0) &&      #  While more words
           ($ARGV[0] !~ m/^-/)) {      #  and not an option word,
      $sampleList{ shift @ARGV } = 1;  #  add sample to the list
    }
  }

  #  --read-type is available for things that operate on reads directly:
  #  fetch, list, read-stats, read-filter
  #
  elsif ((($mode eq "fetch") ||
          ($mode eq "list") ||
          ($mode eq "read-stats") ||
          ($mode eq "read-filter")) && ($arg =~ "--read-types?")) {
    while ((scalar(@ARGV) > 0) &&     #  While more words
           ($ARGV[0] !~ m/^-/)) {     #  and not an option word,
      $readTypes{ shift @ARGV } = 1;  #  add data-type to the list
    }
  }

  #  --submit is available for things that compute.
  #
  elsif ((($mode eq "read-stats") ||
          ($mode eq "read-filter") ||
          ($mode eq "genomescope") ||
          ($mode eq "hapmers") ||
          ($mode eq "yakmers") ||
          ($mode eq "ribotin") ||
          ($mode eq "assemble") ||
          ($mode eq "analyze")) && ($arg eq "--submit")) {
    $opts{"submit"} = 1;
  }


  #  HELP          (no options, though --samples is accepted)
  #  LIST          --samples ...  --files --summary --details --read-type     #  displays list of samples by defailt
  #  STATUS        --samples ...
  #  FETCH         --samples ...
  #  READ-STATS    --samples ...  --submit              --read-type hifi|ont
  #  READ-FILTER   --samples ...  --submit  --cleanup
  #  GENOMESCOPE   --samples ...  --submit
  #  HAPMERS       --samples ...  --submit
  #  YAKMERS       --samples ...  --submit
  #  ASSEMBLE      --samples ...  --submit  --cleanup --archive   [--canu-trio --canu-hifi --verkko-base --verkko-trio --verkko-hi-c --verkko-thic]
  #  ANALYZE       --samples ...
  #
  #  FETCH and READ-STATS run locally.

  #lsif (($mode eq "list") && ())                       {                        }   #  default is to show a list of sample names
  elsif (($mode eq "list") && ($arg eq "--sizes"))      { $opts{"sizes"}    = 1; }   #  total size of aws/downloaded files for each sample
  elsif (($mode eq "list") && ($arg eq "--files"))      { $opts{"files"}    = 1; }   #  same, but also show each file
  elsif (($mode eq "list") && ($arg eq "--summary"))    { $opts{"summary"}  = 1; }   #  total coverage per type and heterozygosity from genomescope
  elsif (($mode eq "list") && ($arg eq "--details"))    { $opts{"details"}  = 1; }   #  coverage per file

  elsif (($mode eq "status") && ($arg eq "--jobs"))     { $opts{"jobs"}     = 1; }
  elsif (($mode eq "status") && ($arg eq "--times"))    { $opts{"times"}    = 1; }

  elsif (($mode eq "status") && ($arg eq "--finished")) { $opts{"finished"} = 1; }
  elsif (($mode eq "status") && ($arg eq "--running"))  { $opts{"running"}  = 1; }
  elsif (($mode eq "status") && ($arg eq "--failed"))   { $opts{"failed"}   = 1; }
  elsif (($mode eq "status") && ($arg eq "--waiting"))  { $opts{"waiting"}  = 1; }

  #lsif (($mode eq "read-stats") && ($arg eq "--submit"))  { $opts{"submit"}  = 1; }

  elsif ($arg eq "--canu-hifi")     { $opts{"flavor"} = "canu-hifi";   }
  elsif ($arg eq "--canu-trio")     { $opts{"flavor"} = "canu-trio";   }

  elsif ($arg eq "--verkko-base")   { $opts{"flavor"} = "verkko-base"; }   #  Correction, MBG up through untip.  All depend on this.
  elsif ($arg eq "--verkko-trio")   { $opts{"flavor"} = "verkko-trio"; }   #  Continue 'base' into a trio assembly.
  elsif ($arg eq "--verkko-hi-c")   { $opts{"flavor"} = "verkko-hi-c"; }   #  Continue 'base' into a hi-c assembly.
  elsif ($arg eq "--verkko-thic")   { $opts{"flavor"} = "verkko-thic"; }   #  Continue 'base' into a trio assembly with hi-c for rdna.

  elsif ($arg eq "--verkko-full")   { $opts{"flavor"} = "verkko-full"; }   #  Do everything in one job.


  elsif (($mode eq "assemble") && ($arg eq "--cleanup")) {
    $opts{"cleanup"} = 1;
  }

  elsif (($mode eq "assemble") && ($arg eq "--archive")) {
    $opts{"archive"} = 1;
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

    push @errs, "Sample$ss $ms not known.";
  }
}



#
#  Expand fetch options and check for invalid ones.
#
if (scalar(keys %readTypes) == 0) {
  $readTypes{'all'} = 1;
}
if (exists($readTypes{'mati'})) {
  $readTypes{'mat-ilmn'} = 1;   delete $readTypes{'mati'};
}
if (exists($readTypes{'pati'})) {
  $readTypes{'pat-ilmn'} = 1;   delete $readTypes{'pati'};
}
if (exists($readTypes{'all'})) {
  $readTypes{$_} = 1   foreach qw(hifi ont hic ilmn mat-ilmn pat-ilmn);
  delete $readTypes{'all'};
}
foreach my $f (keys %readTypes) {   #  Check for invalid fetch options.
  if (($f ne 'hifi') && ($f ne 'hifi-cutadapt') && ($f ne 'ont') &&
      ($f ne 'hic') && ($f ne 'hic1') && ($f ne 'hic2') &&
      ($f ne 'ilmn') && 
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
    ($mode ne "genomescope") &&
    ($mode ne "hapmers") &&
    ($mode ne "yakmers") &&
    ($mode ne "ribotin") &&
    ($mode ne "assemble") &&
    ($mode ne "analyze")) {
  push @errs, "Invalid mode '$mode'.";
}

if (($opts{'submit'}) && ($ENV{'HOSTNAME'} =~ m/helix/)) {
  push @errs, "Can't --submit on host helix.";
}

if (scalar(keys %samples) == 0) {
  push @errs, "No known samples; you probably need to specify --v3 or --v4 or similar.";
}

if (($mode eq "help") || (scalar(@errs) > 0)) {
  print "usage: $0 mode [options]";
  print "  <see the docs>\n";
  print "\n";
  print "ERROR: $_\n"  foreach (@errs);
  exit(1);
}



if ($mode eq "list") {
  if ((exists($opts{"summary"})) ||
      (exists($opts{"details"})) ||
      (exists($opts{"sizes"}))) {
    foreach my $s (sort keys %sampleList) {
      displaySample($s, \%readTypes, \%opts);
    }
  }
  else {
    foreach my $s (sort keys %sampleList) {  #  We never get UNKNOWN samples here; they
      print "$s\n";                          #  generate an error above.
    }
  }
}

elsif ($mode eq "status") {
  foreach my $s (sort keys %sampleList) {
    if ($opts{"flavor"} eq "") {
      $opts{"flavor"} = "canu-trio";    checkStatus($s, \%opts);
      $opts{"flavor"} = "canu-hifi";    checkStatus($s, \%opts);
      $opts{"flavor"} = "verkko-trio";  checkStatus($s, \%opts);
      $opts{"flavor"} = "verkko-hi-c";  checkStatus($s, \%opts);
      $opts{"flavor"} = "verkko-full";  checkStatus($s, \%opts);
      $opts{"flavor"} = "";
    }
    else {
      checkStatus($s, \%opts);
    }
  }
}

elsif ($mode eq "fetch") {
  foreach my $s (sort keys %sampleList) {
    fetchData($s, \%readTypes, \%opts);
  }
}

elsif ($mode eq "read-stats") {
  foreach my $s (sort keys %sampleList) {
    summarizeReads($s, \%readTypes, \%opts);
  }
}

elsif ($mode eq "read-filter") {
  foreach my $s (sort keys %sampleList) {
    filterReads($s, \%readTypes, ,\%opts);
  }
}

elsif ($mode eq "genomescope") {
  foreach my $s (sort keys %sampleList) {
    computeGenomescope($s, \%opts);
  }
}

elsif ($mode eq "hapmers") {
  foreach my $s (sort keys %sampleList) {
    computeHapmers($s, \%opts);
  }
}

elsif ($mode eq "yakmers") {
  foreach my $s (sort keys %sampleList) {
    computeYakmers($s, \%opts);
  }
}

elsif ($mode eq "assemble") {
  my $cleanup = $opts{'cleanup'};
  my $archive = $opts{'archive'};
  my $compute = !$cleanup && !$archive;

  foreach my $s (sort keys %sampleList) {
    cleanupAssembly($s, \%opts)   if ($cleanup);
    archiveAssembly($s, \%opts)   if ($archive);
    computeAssembly($s, \%opts)   if ($compute);
  }
}

elsif ($mode eq "analyze") {
  my $bc = "$root/hprc-cache/busco";     #  Also hardcoded in analyze.sh

  foreach my $db (qw(primates)) {   #  Also included euarchontoglires for some reason
    if (! -e "$bc") {
      system("mkdir -p $bc");
    }
    if (! -e "$bc/${db}_odb10.2024-01-08.tar.gz") {
      print STDERR "Fetching BUSCO ${db}.\n";
      system("cd $bc ; curl -LRO https://busco-data.ezlab.org/v5/data/lineages/${db}_odb10.2024-01-08.tar.gz 2> ${db}.fetch.err");
    }
    if (! -e "$bc/${db}_odb10/dataset.cfg") {
      print STDERR "Unpacking BUSCO ${db} (this is slow).\n";
      system("cd $bc ; tar -xzf ${db}_odb10.2024-01-08.tar.gz");
    }
  }

  foreach my $s (sort keys %sampleList) {
    if (($opts{"flavor"} eq "") || ($opts{"flavor"} eq "verkko-full")) {
      #opts{"flavor"} = "canu-trio";    startChromosomeAssignment($s, \%opts);   startTelomereAnalysis($s, \%opts);    #  These don't support Canu
      #opts{"flavor"} = "canu-hifi";    startChromosomeAssignment($s, \%opts);   startTelomereAnalysis($s, \%opts);    #  output names.
      $opts{"flavor"} = "verkko-trio";  startChromosomeAssignment($s, \%opts);   startTelomereAnalysis($s, \%opts);
      $opts{"flavor"} = "verkko-hi-c";  startChromosomeAssignment($s, \%opts);   startTelomereAnalysis($s, \%opts);
      $opts{"flavor"} = "verkko-thic";  startChromosomeAssignment($s, \%opts);   startTelomereAnalysis($s, \%opts);
      $opts{"flavor"} = "";
    }
    else {
      startChromosomeAssignment($s, \%opts);
      startTelomereAnalysis    ($s, \%opts);
    }
  }
}

elsif ($mode eq "ribotin") {
  foreach my $s (sort keys %sampleList) {
    if (($opts{"flavor"} eq "") || ($opts{"flavor"} eq "verkko-full")) {
      $opts{"flavor"} = "verkko-hi-c";  startRibotinAnalysis($s, \%opts);
      $opts{"flavor"} = "";
    }
    else {
      startRibotinAnalysis($s, \%opts);
    }
  }
}
