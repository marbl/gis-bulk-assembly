
###############################################################################
#
#  This file is part of a pipeline to (help) run Verkko assemblies on
#  numerous HPRC samples.
#
#  This is a 'United States Government Work', and is released in the public
#  domain.
#
##

package hprc::assemble;
require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(startAssembly checkAssembly);

use strict;
use warnings "all";
no  warnings "uninitialized";

use hprc::samples;

my $root = "/data/walenzbp/hprc";


sub getFiles ($$) {
  my $samp  = shift @_;
  my $type  = shift @_;
  my $files = $samples{$samp}{$type};
  my @flist;
  my @missing;

  foreach my $f (@$files) {
    my $awsf = $f;   #  So we don't accidentally corrupt the file list.
    my $locf = $f;
    my $locn = $f;

    $locf =~ s!/!--!g;
    $locf =~ s!^s3:----human-pangenomics--submissions--!$root/aws-data/$samp/!;

    $locn =~ s!^s3://human-pangenomics/submissions/!!;

    if (-e $locf) {
      #printf STDERR "OK      %8s:%s\n", $type, $locn;
      push @flist, $locf;
    } else {
      push @missing, "$type\0$locn";
    }
  }

  my $flist = join " \\\n         ", @flist;

  if (scalar(@missing) > 0) {
    print STDERR "Can't launch the assembly, inputs missing.\n";

    foreach my $m (@missing) {
      my ($t, $f) = split '\0', $m;

      printf STDERR "  %8s - %s\n", $t, $f;
    }

    exit(1);
  }

  return($flist);
}


sub startAssembly ($$) {
  my $samp   = shift @_;
  my $submit = shift @_;

  #  Make a place to work.

  if (! -d "$root/assemblies/$samp") {
    system("mkdir -p $root/assemblies/$samp");
  }

  #  Check that inputs exist.

  my $hifi = getFiles($samp, "hifi");
  my $nano = getFiles($samp, "ont");

  #  Make sure it isn't running.

  if (-e "$root/assemblies/$samp.jid") {
    my ($jobid, $state, $reason, $out) = (undef, undef, undef, undef);

    open(JID, "< $root/assemblies/$samp.jid") or die "Failed to open '$root/assemblies/$samp.jid' for reading: $!\n";
    my $jid = int(<JID>);
    close(JID);

    open(JOBS, "dashboard_cli jobs --noheader --tab --fields jobid,state,state_reason,std_out --jobid $jid|") or die "Failed to open 'dashboard_cli jobs' for reading: $!\n";
    while (<JOBS>) {
      ($jobid, $state, $reason, $out) = split '\s+', $_;

      last if ($jobid eq $jid);
    }
    close(JOBS);

    if (!defined($jobid)) {
      print STDERR "Didn't get a report for job $jid from Slurm.  Can't check if I can safely submit.\n";
      exit(0);
    }

    if    (($state eq "CONFIGURING") ||
           ($state eq "PENDING")) {
      print STDERR "Job $jobid is in state $state.  Cannot submit it again.\n";
      exit(0);
    }
    elsif (($state eq "RUNNING") ||
           ($state eq "COMPLETING") ||
           ($state eq "STOPPED") ||
           ($state eq "SUSPENDED")) {
      print STDERR "Job $jobid is in state $state.  Cannot submit it again.\n";
      exit(0);
    }
    elsif (($state eq "COMPLETED") ||
           ($state eq "CANCELLED") ||
           ($state eq "FAILED") ||
           ($state eq "TIMEOUT")) {
      #print STDERR "Job $jobid is in state $state.  Will submit again.\n";
    }
    else {
      print STDERR "Job $jobid is in state $state.  Not sure if I can safely submit again.\n";
      exit(0);
    }
  }

  #  Launch.

  open(CMD, "> $root/assemblies/$samp.sh") or die "Failed to open 'assemblies/$samp.sh' for writing: $!\n";
  print CMD "#!/bin/sh\n";
  print CMD "#\n";
  print CMD "#SBATCH --cpus-per-task=2\n";
  print CMD "#SBATCH --mem=16g\n";
  print CMD "#SBATCH --time=4-0\n";
  print CMD "#SBATCH --output=$root/assemblies/$samp.%j.log\n";
  print CMD "#SBATCH --job-name=verkko$samp\n";
  print CMD "#\n";
  print CMD "cd $root/assemblies/$samp\n";
  print CMD "\n";
  print CMD "module load winnowmap\n";
  print CMD "module load mashmap\n";
  print CMD "module load samtools\n";
  print CMD "\n";
  print CMD "echo 'PYTHON VERSION:'\n";
  print CMD "command -v python\n";
  print CMD "python --version\n";
  print CMD "\n";
  print CMD "echo ''\n";
  print CMD "echo 'PYTHON LIBRARIES:'\n";
  print CMD "echo  > x.py 'from importlib import metadata'\n";
  print CMD "echo >> x.py 'for dist in metadata.distributions():'\n";
  print CMD "echo >> x.py '  print(dist._path)'\n";
  print CMD "python x.py\n";
  print CMD "\n";
  print CMD "echo ''\n";
  print CMD "echo 'NETWORKX VERSION'\n";
  print CMD "echo  > x.py 'from importlib.metadata import version'\n";
  print CMD "echo >> x.py 'version(\\\"networkx\\\")'\n";
  print CMD "python x.py\n";
  print CMD "\n";
  print CMD "echo ''\n";
  print CMD "echo 'NETWORKX TEST'\n";
  print CMD "echo  > x.py 'import networkx as nx'\n";
  print CMD "echo >> x.py 'G = nx.Graph()'\n";
  print CMD "echo >> x.py 'print(G)'\n";
  print CMD "python x.py\n";
  print CMD "\n";
  print CMD "rm x.py\n";
  print CMD "\n";
  print CMD "$root/software/verkko/bin/verkko --slurm -d . \\\n";
  print CMD "  --ovb-run 8 32 8 \\\n";
  print CMD "  --hifi $hifi \\\n";
  print CMD "  --nano $nano\n";
  print CMD "\n";
  close(CMD);

  print STDOUT "sbatch $root/assemblies/$samp.sh > $root/assemblies/$samp.jid\n";
  system("sbatch $root/assemblies/$samp.sh > $root/assemblies/$samp.jid 2>&1")   if ($submit);
}


sub checkAssembly ($$) {
}

1;
