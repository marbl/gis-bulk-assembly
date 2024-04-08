
###############################################################################
#
#  This file is part of a pipeline to (help) run Verkko assemblies on
#  numerous HPRC samples.
#
#  This is a 'United States Government Work', and is released in the public
#  domain.
#
##

package hprc::genomescope;
require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(computeGenomescope);

use strict;
use warnings "all";
no  warnings "uninitialized";

use hprc::samples;
use hprc::assemble;

my $genomescope = "/data/Phillippy/tools/genomescope2.0/genomescope.R";

sub computeGenomescope ($$) {
  my $samp = shift @_;
  my $opts = shift @_;
  my $geno   = "$root/aws-data/$samp/genomescope";

  #  Check if outputs exist.

  if ((-e "$geno/hifi/summary.txt") &&
      (-e "$geno/ilmn/summary.txt") &&
      (-e "$geno/mati/summary.txt") &&
      (-e "$geno/pati/summary.txt")) {
    print "Outputs $samp/hifi/ and\n";
    print "        $samp/ilmn/ and\n";
    print "        $samp/mati/ and\n";
    print "        $samp/pati/ exist, nothing to do.\n";
    return;
  }

  #  Fail early if genomescope isn't available.

  if (! -e $genomescope) {
    print STDERR "ERROR: genomescope not found at:\n";
    print STDERR "           '$genomescope'\n";
    exit(1);
  }

  #  Make a place to work.

  if (! -d "$geno") {
    system("mkdir -p $geno");
  }

  #  Check that inputs exist.

  my $hifi = getFiles($samp, "hifi-cutadapt");
  my $ilmn = getFiles($samp, "ilmn");
  my $mati = getFiles($samp, "mat-ilmn");
  my $pati = getFiles($samp, "pat-ilmn");

  #  Make sure it isn't running.

  if (-e "$geno/genomescope.jid") {
    print "$samp is currently running.\n";
    return;
  }

  my $meryl = "$root/software/meryl/build/bin/meryl";
  my $d     = "/lscratch/\$SLURM_JOB_ID";

  #  Launch.
  #    HG04199 used just over 250GB disk.

  open(CMD, "> $geno/genomescope.sh") or die "Failed to open '$geno/genomescope.sh' for writing: $!\n";
  print CMD "#!/bin/sh\n";
  print CMD "#\n";
  print CMD "#SBATCH --gres=lscratch:512\n";
  print CMD "#SBATCH --cpus-per-task=16\n";
  print CMD "#SBATCH --mem=128g\n";
  print CMD "#SBATCH --time=1-0\n";
  print CMD "#SBATCH --output=$geno/genomescope.err\n";
  print CMD "#SBATCH --job-name=gscope$samp\n";
  print CMD "#\n";
  print CMD "set -x\n";
  print CMD "\n";
  print CMD "export REF_CACHE=$root/samtools-ref-cache\n";
  print CMD "\n";
  print CMD "module load R\n";
  print CMD "\n";
  print CMD "cd $geno\n";
  print CMD "\n";
  print CMD "if [ ! -e hifi.k21.histogram ] ; then\n";
  print CMD "  $meryl k=21 threads=16 memory=96g count \\\n";
  print CMD "    output $d/hifi.meryl \\\n";
  print CMD "         $hifi \\\n";
  print CMD "  && \\\n";
  print CMD "  $meryl histogram $d/hifi.meryl > hifi.k21.histogram\n";
  print CMD "  rm -rf $d/hifi.meryl\n";
  print CMD "fi\n";
  print CMD "\n";
  print CMD "if [ ! -e ilmn.k21.histogram ] ; then\n";
  print CMD "  $meryl k=21 threads=16 memory=96g count \\\n";
  print CMD "    output $d/ilmn.meryl \\\n";
  print CMD "         $ilmn \\\n";
  print CMD "  && \\\n";
  print CMD "  $meryl histogram $d/ilmn.meryl > ilmn.k21.histogram\n";
  print CMD "  rm -rf $d/ilmn.meryl\n";
  print CMD "fi\n";
  print CMD "\n";
  print CMD "if [ ! -e mati.k21.histogram ] ; then\n";
  print CMD "  $meryl k=21 threads=16 memory=96g count \\\n";
  print CMD "    output $d/mati.meryl \\\n";
  print CMD "         $mati \\\n";
  print CMD "  && \\\n";
  print CMD "  $meryl histogram $d/mati.meryl > mati.k21.histogram\n";
  print CMD "  rm -rf $d/mati.meryl\n";
  print CMD "fi\n";
  print CMD "\n";
  print CMD "if [ ! -e pati.k21.histogram ] ; then\n";
  print CMD "  $meryl k=21 threads=16 memory=96g count \\\n";
  print CMD "    output $d/pati.meryl \\\n";
  print CMD "         $pati \\\n";
  print CMD "  && \\\n";
  print CMD "  $meryl histogram $d/pati.meryl > pati.k21.histogram\n";
  print CMD "  rm -rf $d/pati.meryl\n";
  print CMD "fi\n";
  print CMD "\n";
  print CMD "\n";
  print CMD "if [ -e hifi.k21.histogram ] ; then\n";
  print CMD "  $genomescope -i hifi.k21.histogram -o hifi -p 2 -k 21\n";
  print CMD "fi\n";
  print CMD "\n";
  print CMD "if [ -e hifi.k21.histogram ] ; then\n";
  print CMD "  $genomescope -i ilmn.k21.histogram -o ilmn -p 2 -k 21\n";
  print CMD "fi\n";
  print CMD "\n";
  print CMD "if [ -e hifi.k21.histogram ] ; then\n";
  print CMD "  $genomescope -i mati.k21.histogram -o mati -p 2 -k 21\n";
  print CMD "fi\n";
  print CMD "\n";
  print CMD "if [ -e hifi.k21.histogram ] ; then\n";
  print CMD "  $genomescope -i pati.k21.histogram -o pati -p 2 -k 21\n";
  print CMD "fi\n";
  print CMD "\n";
  print CMD "rm -f $geno/genomescope.jid\n";
  print CMD "\n";
  print CMD "exit 0\n";
  close(CMD);

  print STDOUT "sbatch $geno/genomescope.sh > $geno/genomescope.jid\n";
  system("sbatch $geno/genomescope.sh > $geno/genomescope.jid 2>&1")   if (exists($$opts{"submit"}));
}

1;
