
###############################################################################
#
#  This file is part of a pipeline to (help) run Verkko assemblies on
#  numerous HPRC samples.
#
#  This is a 'United States Government Work', and is released in the public
#  domain.
#
##

package hprc::yakmers;
require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(computeYakmers);

use strict;
use warnings "all";
no  warnings "uninitialized";

use hprc::samples;
use hprc::assemble;


sub emitYak ($$$$) {
  my $yako = shift @_;
  my $samp = shift @_;
  my $type = shift @_;
  my $inps = shift @_;

  my $p1 = "$yako/yak-$samp-$type-1.fifo";
  my $p2 = "$yako/yak-$samp-$type-2.fifo";

  print CMD "\n";
  print CMD "\n";
  print CMD "if [ ! -e $yako/$type.yak ] ; then\n";
  print CMD "  mkfifo $p1\n";
  print CMD "  mkfifo $p2\n";
  print CMD "\n";
  print CMD "  $root/software/seqrequester/build/bin/seqrequester extract -fasta $inps > $p1 &\n";
  print CMD "  $root/software/seqrequester/build/bin/seqrequester extract -fasta $inps > $p2 &\n";
  print CMD "\n";
  print CMD "  yak count -k31 -b37 -t \$yakCPUs -o $yako/$type.yak $p1 $p2 \\\n";
  print CMD "  || \\\n";
  print CMD "  rm -f $yako/$type.yak\n";
  print CMD "\n";
  print CMD "  rm -f $p1\n";
  print CMD "  rm -f $p2\n";
  print CMD "else\n";
  print CMD "  echo $yako/$type.yak exists, NOT recomputing.\n";
  print CMD "fi\n";
}


sub computeYakmers ($$) {
  my $samp   = shift @_;
  my $opts   = shift @_;
  my $yako   = "$root/aws-data/$samp/yak";   #  Yakko waves to Wakko and Dot.

  print "\n";
  print "$samp\n";

  #  Check if outputs exist.

  system("mkdir -p $yako");

  if ((-e "$yako/ilmn.yak") &&
      (-e "$yako/mati.yak") &&
      (-e "$yako/pati.yak")) {
    print "Outputs $yako/ilmn.yak and\n";
    print "        $yako/mati.yak and\n";
    print "        $yako/pati.yak exist, nothing to do.\n";
    return;
  }

  #  Check that inputs exist.

  if ((numFiles($samp, "ilmn")     == 0) ||
      (numFiles($samp, "mat-ilmn") == 0) ||
      (numFiles($samp, "pat-ilmn") == 0)) {
    print "  No child Illumina data exists.\n"      if (numFiles($samp, "ilmn") == 0);
    print "  No maternal Illumina data exists.\n"   if (numFiles($samp, "mat-ilmn") == 0);
    print "  No paternal Illumina data exists.\n"   if (numFiles($samp, "pat-ilmn") == 0);
    return;
  }

  my $ilmn = getFiles($samp, "ilmn");
  my $mati = getFiles($samp, "mat-ilmn");
  my $pati = getFiles($samp, "pat-ilmn");

  if (($ilmn eq "") ||
      ($mati eq "") ||
      ($pati eq "")) {
    print STDERR "  Inputs missing (see above; check 'hprc.pl list --sample $samp').\n";
    return;
  }

  #  Make sure it isn't running.

  if (-e "$yako/yak.jid") {
    print "$samp is currently running.\n";
    return;
  }
  
  #  Launch.

  open(CMD, "> $yako/yak.sh") or die "Failed to open '$yako/yak.sh' for writing: $!\n";
  print CMD "#!/bin/sh\n";
  print CMD "#\n";
  #rint CMD "#SBATCH --gres=lscratch:512\n";
  print CMD "#SBATCH --cpus-per-task=16\n";
  print CMD "#SBATCH --mem=96g\n";
  print CMD "#SBATCH --time=8:00:00\n";
  print CMD "#SBATCH --output=$yako/yak.err\n";
  print CMD "#SBATCH --job-name=yak$samp\n";
  print CMD "#\n";
  print CMD "set -x\n";
  print CMD "\n";
  print CMD "export REF_CACHE=$root/samtools-ref-cache\n";
  print CMD "\n";
  print CMD "module load yak\n";
  print CMD "\n";
  print CMD "cd $yako\n";
  print CMD "\n";
  print CMD "yak version\n";
  print CMD "\n";

  print CMD "yakCPUs=\$SLURM_JOB_CPUS_PER_NODE\n";
  print CMD "\n";
  print CMD "if [ x\$yakCPUs = x ] ; then\n";
  print CMD "  echo 'WARNING: SLURM_JOB_CPUS_PER_NODE not set, using minimal CPUs instead.'\n";
  print CMD "  yakCPUs=2\n";
  print CMD "fi\n";

  emitYak($yako, $samp, "ilmn", $ilmn);
  emitYak($yako, $samp, "mati", $mati);
  emitYak($yako, $samp, "pati", $pati);

  print CMD "\n";
  print CMD "\n";
  print CMD "unlink $yako/yak.jid\n";
  print CMD "\n";
  print CMD "exit 0\n";
  close(CMD);

  print STDOUT "sbatch $yako/yak.sh > $yako/yak.jid\n";
  system("sbatch $yako/yak.sh > $yako/yak.jid 2>&1")   if (exists($$opts{"submit"}));
}

1;
