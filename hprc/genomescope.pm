
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

use hprc::aws;
use hprc::samples;
use hprc::assemble;

my $genomescope = "/data/Phillippy/tools/genomescope2.0/genomescope.R";
my $meryl       = "$rsoft/meryl/build/bin/meryl";

sub computeGenomescope ($$) {
  my $samp  = shift @_;
  my $opts  = shift @_;
  my $geno  = "$data/$samp/genomescope";
  my $locd  = "/lscratch/\$SLURM_JOB_ID";

  #  Fail early if genomescope isn't available.

  if (! -e $genomescope) {
    print STDERR "ERROR: genomescope not found at:\n";
    print STDERR "           '$genomescope'\n";
    exit(1);
  }

  #  Check if outputs exist.

  my $compl = ((-e "$geno/hifi/summary.txt") &&
               (-e "$geno/ilmn/summary.txt") &&
               (-e "$geno/mati/summary.txt") &&
               (-e "$geno/pati/summary.txt"));

  #  Check that inputs exist.

  my ($nhifi, $hifi) = (numFiles($samp, "hifi-cutadapt"), getDownloadedFiles($samp, "hifi-cutadapt"));
  my ($nilmn, $ilmn) = (numFiles($samp, "ilmn"),          getDownloadedFiles($samp, "ilmn"));
  my ($nmati, $mati) = (numFiles($samp, "mat-ilmn"),      getDownloadedFiles($samp, "mat-ilmn"));
  my ($npati, $pati) = (numFiles($samp, "pat-ilmn"),      getDownloadedFiles($samp, "pat-ilmn"));

  my $downl = ((($nhifi == 0) || ($hifi ne "")) &&   #  True if there are no reads in the data set or
               (($nilmn == 0) || ($ilmn ne "")) &&   #  all reads are downloaded.  getDownloadedFiles()
               (($nmati == 0) || ($mati ne "")) &&   #  returns empty-string if any file is missing.
               (($npati == 0) || ($pati ne "")));

  #  Create a script - if data is donwloaded and script doesn't exist.

  if ($downl && ! -e "$geno/genomescope.sh") {
    system("mkdir -p $geno")   if (! -d "$geno");

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
    print CMD "export REF_CACHE=$ENV{'REF_CACHE'}\n";
    print CMD "\n";
    print CMD "module load R\n";
    print CMD "\n";
    print CMD "cd $geno\n";
    print CMD "\n";
    print CMD "if [ ! -e hifi.k21.histogram ] ; then\n";
    print CMD "  $meryl k=21 threads=16 memory=96g count \\\n";
    print CMD "    output $locd/hifi.meryl \\\n";
    print CMD "         $hifi \\\n";
    print CMD "  && \\\n";
    print CMD "  $meryl histogram $locd/hifi.meryl > hifi.k21.histogram\n";
    print CMD "  rm -rf $locd/hifi.meryl\n";
    print CMD "fi\n";
    print CMD "\n";
    print CMD "if [ ! -e ilmn.k21.histogram ] ; then\n";
    print CMD "  $meryl k=21 threads=16 memory=96g count \\\n";
    print CMD "    output $locd/ilmn.meryl \\\n";
    print CMD "         $ilmn \\\n";
    print CMD "  && \\\n";
    print CMD "  $meryl histogram $locd/ilmn.meryl > ilmn.k21.histogram\n";
    print CMD "  rm -rf $locd/ilmn.meryl\n";
    print CMD "fi\n";
    print CMD "\n";
    print CMD "if [ ! -e mati.k21.histogram ] ; then\n";
    print CMD "  $meryl k=21 threads=16 memory=96g count \\\n";
    print CMD "    output $locd/mati.meryl \\\n";
    print CMD "         $mati \\\n";
    print CMD "  && \\\n";
    print CMD "  $meryl histogram $locd/mati.meryl > mati.k21.histogram\n";
    print CMD "  rm -rf $locd/mati.meryl\n";
    print CMD "fi\n";
    print CMD "\n";
    print CMD "if [ ! -e pati.k21.histogram ] ; then\n";
    print CMD "  $meryl k=21 threads=16 memory=96g count \\\n";
    print CMD "    output $locd/pati.meryl \\\n";
    print CMD "         $pati \\\n";
    print CMD "  && \\\n";
    print CMD "  $meryl histogram $locd/pati.meryl > pati.k21.histogram\n";
    print CMD "  rm -rf $locd/pati.meryl\n";
    print CMD "fi\n";
    print CMD "\n";
    print CMD "\n";
    print CMD "if [ -e hifi.k21.histogram ] ; then\n";
    print CMD "  rm -f hifi hifi.FAIL\n";
    print CMD "  $genomescope -i hifi.k21.histogram -o hifi -p 2 -k 21\\\n";
    print CMD "  ||\\\n";
    print CMD "  mv -f hifi hifi.FAIL\n";
    print CMD "fi\n";
    print CMD "\n";
    print CMD "if [ -e ilmn.k21.histogram ] ; then\n";
    print CMD "  rm -f ilmn ilmn.FAIL\n";
    print CMD "  $genomescope -i ilmn.k21.histogram -o ilmn -p 2 -k 21\n";
    print CMD "  ||\\\n";
    print CMD "  mv -f ilmn ilmn.FAIL\n";
    print CMD "fi\n";
    print CMD "\n";
    print CMD "if [ -e mati.k21.histogram ] ; then\n";
    print CMD "  rm -f mati mati.FAIL\n";
    print CMD "  $genomescope -i mati.k21.histogram -o mati -p 2 -k 21\n";
    print CMD "  ||\\\n";
    print CMD "  mv -f mati mati.FAIL\n";
    print CMD "fi\n";
    print CMD "\n";
    print CMD "if [ -e pati.k21.histogram ] ; then\n";
    print CMD "  rm -f pati pati.FAIL\n";
    print CMD "  $genomescope -i pati.k21.histogram -o pati -p 2 -k 21\n";
    print CMD "  ||\\\n";
    print CMD "  mv -f pati pati.FAIL\n";
    print CMD "fi\n";
    print CMD "\n";
    print CMD "rm -f $geno/genomescope.jid\n";
    print CMD "\n";
    print CMD "exit 0\n";
    close(CMD);
  }

  if    ($compl)                          { print "$samp/genomescope - FINISHED\n"; }
  elsif (-e "$geno/genomescope.jid")      { print "$samp/genomescope - RUNNING\n"; }
  elsif (-e "$geno/genomescope.err")      { print "$samp/genomescope - CRASHED\n"; }
  elsif (! $downl)                        { print "$samp/genomescope - NOT-FETCHED\n"; }
  elsif (! $$opts{"submit"})              { print "$samp/genomescope - READY-TO-COMPUTE\n"; }
  else                                    { print "$samp/genomescope - SUBMITTED\n"; system("sbatch $geno/genomescope.sh > $geno/genomescope.jid"); }
}

1;
