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
#use hprc::status;
use hprc::assemble;


my $mode   = undef;
my $sample = undef;
my %fetch;

my $help   = undef;
my @errs;

while (scalar(@ARGV) > 0) {
    my $arg = $ARGV[0];  shift @ARGV;

    if    ($arg eq "--help")    { $help = "help";          }
    elsif ($arg eq "--list")    { $mode = "list";          }
    elsif ($arg eq "--status")  { $mode = "status";        }
    elsif ($arg eq "--sample")  { $sample = shift @ARGV;   }

    elsif ($arg eq "--fetch")   {
        $mode  = "fetch";

        if ((scalar(@ARGV) == 0) ||     #  If only '--fetch' and no
            ($ARGV[0] =~ m/^-/)) {      #  data-type, default to 'all'.
            $fetch{'all'} = 1;
        }

        while ((scalar(@ARGV) > 0) &&   #  While more words
               ($ARGV[0] !~ m/^-/)) {   #  and not an option word,
            $fetch{ shift @ARGV } = 1;  #  add data-type to list of fetches.
        }

        if (exists($fetch{'all'})) {
            $fetch{$_} = 1   foreach qw(hifi ont hic ilmn mat-ilmn pat-ilmn);
            delete $fetch{'all'};
        }
    }

    elsif ($arg eq "--read-stats") {
    }

    elsif ($arg eq "--trio") {
    }

    elsif ($arg eq "--assemble") {
      $mode = "assemble";
    }

    elsif ($arg eq "--post") {
    }

    else {
        push @errs, "unknown option '$arg'";
        $help = 1;
    }
}

foreach my $f (keys %fetch) {   #  Check for invalid fetch options.
  if (($f ne 'hifi') && ($f ne 'ont') && ($f ne 'hic') && ($f ne 'ilmn') &&
      ($f ne 'mat-ilmn') && ($f ne 'pat-ilmn')) {
    delete $fetch{$_}  foreach qw(hifi ont hic ilmn mat-ilmn pat-ilmn);
    push @errs, "Invalid '--fetch' options '" . join("', '", sort keys %fetch) . "'.";
    push @errs, "  Must be one of:";
    push @errs, "    'all'                         (everything)";
    push @errs, "    'hifi', 'ont', 'hic', 'ilmn', (child)";
    push @errs, "    'mat-ilmn', 'pat-ilmn'        (parental)";
    last;
  }
}

if (($help ne undef) || ($mode eq undef) || (scalar(@errs) > 0)) {
    print "usage: $0 [...]";
    print "  --list              - list all known samples\n";
    print "  --status            - show status of all samples\n";
    print "  --sample [sample]   - restrict operation to a single sample\n";
    print "\n";
    print "ERROR: $_\n"  foreach (@errs);
    exit(1);
}


#oadSamples("batch-1-test.tsv");
loadSamples("batch-2-with-hic.tsv");

if (defined($sample) && !exists($samples{$sample})) {
    die "ERROR: sample '$sample' not known.\n";
}




if ($mode eq "list") {
  foreach my $s (($sample)  or  (sort keys %samples)) {
    displaySample($s);
  }
}

elsif ($mode eq "status") {
}

elsif ($mode eq "fetch") {
  foreach my $s (($sample)  or  (sort keys %samples)) {
    foreach my $f (sort keys %fetch) {
      fetchData($s, $f);
    }
  }
}

elsif ($mode eq "read-stats") {
}

elsif ($mode eq "status") {
}

elsif ($mode eq "assemble") {
  foreach my $s (($sample)  or  (sort keys %samples)) {
    startAssembly($s);
  }
}

else {
}
