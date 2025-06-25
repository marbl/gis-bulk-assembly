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
use hprc::readCorrect;
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

my $config_file = "$root/config.ini";
my %config;
my $current_section;
open my $fh, '<', $config_file or die "Cannot open config file: $!";
while (<$fh>) {
    chomp;
    next if /^\s*$/ || /^\s*#/;  # skip empty lines and comments

    if (/^\[(\w+)\]/) {
        $current_section = $1;
    }
    elsif ($current_section && /^(\w+)\s*=\s*(.+)$/) {
        my ($key, $value) = ($1, $2);
        # Remove quotes if present
        $value =~ s/^['"]|['"]$//g;
        # Parse list if it's a list (e.g., ['a','b'])
        if ($value =~ /^\[.*\]$/) {
            my @list = $value =~ /['"]([^'"]+)['"]/g;  # capture the content inside single or double quotes
            @list = map {
              my $item = $_;
              $item =~ s/\$root/$root/g if defined $root;
              $item;
            } @list;
            $config{$current_section}{$key} = \@list;
        } else {
            $value =~ s/\$root/$root/g if defined $root;
            $config{$current_section}{$key} = $value;
        }
    }
}
close $fh;

while (scalar(@ARGV) > 0) {
  my $arg = $ARGV[0];  shift @ARGV;

  #  Select a set of samples to use.
  #
  if ($arg =~ /^--v(\d+)$/) {
    my $v = $1;

   if (exists $config{$v}) {
      my %vars = %{ $config{$v}};
      print STDERR "Configuration: v$v\n";
      foreach my $sample (@{$vars{samples} }) {
         $sample =~ s/'//g;  # remove all single quotes

         if ( -e $sample) {
            loadSamples($sample);
         } else {
            die "Invalid sample file $sample provided, check $config_file\n";
         }
         print STDERR "Sample file:  $sample\n";
      }
      $rasm  = makeAbsolute($vars{rasm}) or die "No assembly output defined for $v\n";
      $rsoft = makeAbsolute($vars{rsoft}) or die "No software defined for $v\n";
      $data  = makeAbsolute($vars{data}) or die "No data location defined for $v\n";
      $refn  = makeAbsolute($vars{refn}) or die "No reference defined for $v\n";
      $refc  = makeAbsolute($vars{refc}) or die "No reference (HPC) defined for $v\n";
      $odb   = $vars{odb} or die " No ODB database defined for $v\n";

      # make sure the paths exist
      if (! -d $root)  { die "Invalid root : $root, check your config\n"; }
      if (! -d $rsoft) { die "Invalid rsoft: $rsoft, check your config\n"; }
      if (! -e $refn)  { die "Invalid ref  : $refn, check your config\n"; }
      if (! -e $refc)  { die "Invalid refc : $refc, check your config\n"; }

      $ENV{'HPRC_ROOT_REFERENCE'}     = "$refn";
      $ENV{'HPRC_ROOT_REFERENCE_HPC'} = "$refc";
      $ENV{'HPRC_ROOT_REFERENCE_ODB'} = "$odb";
      $ENV{'HPRC_ROOT_DATA'}          = "$data";

      print STDERR "Data:         $data\n";
      print STDERR "Software:     $rsoft\n";
      print STDERR "Output asm:   $rasm\n";
      print STDERR "Ref:          $refn\n";
      print STDERR "Ref (HPC):    $refc\n";
      print STDERR "ODB DB:       $odb\n";
    } else {
      die "Unknown configuration v$v requested, check your config file $config_file\n";
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
          ($mode eq "read-correct") ||
          ($mode eq "genomescope") ||
          ($mode eq "hapmers") ||
          ($mode eq "yakmers") ||
          ($mode eq "ribotin") ||
          ($mode eq "assemble") ||
          ($mode eq "analyze")) && ($arg eq "--submit")) {
    $opts{"submit"} = 1;
  }


  #  HELP          (no options, though --samples is accepted)
  #  LIST          --samples ...  --sizes [--files [--s3-paths]] --summary --details --read-type     #  displays list of samples by default
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
  elsif (($mode eq "list") && ($arg eq "--files"))      { $opts{"files"}    = 1; }   #  same, but also show each file (must still supply --sizes)
  elsif (($mode eq "list") && ($arg eq "--s3-paths"))   { $opts{"s3-paths"} = 1; }   #  same, but show the aws s3 path instead of the local path (must still supply --files)
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

  elsif ($arg eq "--hifiasm-trio")  { $opts{"flavor"} = "hifiasm-trio";}
  elsif ($arg eq "--hifiasm-hi-c")  { $opts{"flavor"} = "hifiasm-hi-c";}

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
  if (($f ne 'hifi') && ($f ne 'hifi-cutadapt') && ($f ne 'ont') && ($f ne 'ont-r10') &&
      ($f ne 'hic') && ($f ne 'hic1') && ($f ne 'hic2') &&
      ($f ne 'ilmn') &
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
    ($mode ne "read-correct") &&
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
  push @errs, "No known samples; you probably need to specify --vX based on what you specified in the config.ini file.";
}

if (($mode eq "help") || (scalar(@errs) > 0)) {
  print "usage: $0 mode [options]";
  if ($mode eq "help") {
     open my $fh, '<', "$root/hprc.txt" or die "Cannot open file '$root/hprc.txt': $!";
     print while <$fh>;
     close $fh;
  }
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
    filterReads($s, \%readTypes,\%opts);
  }
}

elsif ($mode eq "read-correct") {
  foreach my $s (sort keys %sampleList) {
    correctReads($s, \%opts);
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
