
###############################################################################
#
#  This file is part of a pipeline to (help) run Verkko assemblies on
#  numerous HPRC samples.
#
#  This is a 'United States Government Work', and is released in the public
#  domain.
#
##

package hprc::analyze;
require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(startChromosomeAssignment startTelomereAnalysis);

use strict;
use warnings "all";
no  warnings "uninitialized";

use hprc::samples;
use hprc::assemble;


sub startChromosomeAssignment ($$) {
  my $samp   = shift @_;
  my $opts   = shift @_;
  my $flav   = $$opts{"flavor"};
  my $submit = $$opts{"submit"};
  my $diro   = "$rasm/$samp/$flav/chr-assign";

  #  Check if outputs exist.

  my $finished = -e "$diro/chr-assign.out";

  #  Check that inputs exist.

  my $ready   = isFinished($samp, $flav);
  my $unavail = dataAvailable($samp, $flav);

  #  Create a script (we don't need to since it exists elsewhere).

  if ($ready && !$finished) {
    system("mkdir -p $diro")   if (! -d "$diro");

    open(CMD, "> $diro/chr-assign.sh") or die "Failed to open '$diro/chr-assign.sh' for writing: $!\n";
    print CMD "#!/bin/sh\n";
    print CMD "#\n";
    print CMD "#SBATCH --cpus-per-task=4\n";
    print CMD "#SBATCH --partition=norm,quick\n";
    print CMD "#SBATCH --mem=4g\n";
    print CMD "#SBATCH --time=1:00:00\n";       #  Actually like 5 minutes
    print CMD "#SBATCH --output=./chr-assign.err\n";
    print CMD "#SBATCH --job-name=va$samp\n";
    print CMD "#\n";
    print CMD "set -e\n";
    print CMD "set -x\n";
    print CMD "\n";
    print CMD "mkdir -p $diro\n";
    print CMD "cd       $diro\n";
    print CMD "\n";
    print CMD "export REF_CACHE=$ENV{'REF_CACHE'}\n";
    print CMD "\n";
    print CMD "module load mashmap\n";
    print CMD "\n";
    print CMD "if [ ! -e asm-chr.mashmap.out ] ; then\n";
    print CMD "  mashmap \\\n";
    print CMD "    --ref /data/Phillippy/references/T2T-CHM13/chm13v2.0.fa \\\n";
    print CMD "    --query ../assembly.fasta \\\n";
    print CMD "    --perc_identity 95 \\\n";
    print CMD "    --segLength 100000 \\\n";
    print CMD "    --threads 4 \\\n";
    print CMD "    --output asm-chr.mashmap.out\n";
    print CMD "fi\n";
    print CMD "\n";
    print CMD "perl $root/hprc/chromosome-assign.pl \\\n";
    print CMD " < asm-chr.mashmap.out \\\n";
    print CMD " > chr-assign.out \\\n";
    print CMD "|| \\\n";
    print CMD "mv chr-assign.out.FAIL\n";
    print CMD "\n";
    print CMD "rm -f chr-assign.jid\n";
    print CMD "\n";
    print CMD "exit 0\n";
  }

  #  Run, if not finished already.

  if    ($finished)                   { print "$samp/$flav/chr-assign - FINISHED\n"; }
  elsif (-e "$diro/chr-assign.jid")   { print "$samp/$flav/chr-assign - RUNNING\n"; }
  elsif (-e "$diro/chr-assign.err")   { print "$samp/$flav/chr-assign - CRASHED\n"; }
  elsif ($unavail)                    { print "$samp/$flav/chr-assign - UNAVAILABLE ($unavail)\n"; }
  elsif (! $ready)                    { print "$samp/$flav/chr-assign - ASSEMBLY-NOT-READY\n"; }
  elsif (! $$opts{"submit"})          { print "$samp/$flav/chr-assign - READY-TO-COMPUTE\n"; }
  else                                { print "$samp/$flav/chr-assign - SUBMITTED\n"; system("sbatch $diro/chr-assign.sh > $diro/chr-assign.jid"); }
}


sub startTelomereAnalysis ($$) {
  my $samp   = shift @_;
  my $opts   = shift @_;
  my $flav   = $$opts{"flavor"};
  my $submit = $$opts{"submit"};
  my $diro   = "$rasm/$samp/$flav/analysis";

  #  Check if outputs exist.

  my $finished = -e "$diro/analysis.complete";

  #  Check that inputs exist.

  my $ready   = isFinished($samp, $flav);
  my $unavail = dataAvailable($samp, $flav);

  #  Create a script (we don't need to since it exists elsewhere).

  if ($ready && !$finished) {
    system("mkdir -p $diro")   if (! -d "$diro");
  }

  #  Run, if not finished already.

  if    ($finished)                   { print "$samp/$flav/analysis   - FINISHED\n"; }
  elsif (-e "$diro/analysis.jid")     { print "$samp/$flav/analysis   - RUNNING\n"; }
  elsif (-e "$diro/analysis.err")     { print "$samp/$flav/analysis   - CRASHED (move $diro/analysis.err to restart)\n"; }
  elsif ($unavail)                    { print "$samp/$flav/analysis   - UNAVAILABLE ($unavail)\n"; }
  elsif (! $ready)                    { print "$samp/$flav/analysis   - ASSEMBLY-NOT-READY\n"; }
  elsif (! $$opts{"submit"})          { print "$samp/$flav/analysis   - READY-TO-COMPUTE\n"; }
  else                                { print "$samp/$flav/analysis   - SUBMITTED\n"; system("sbatch -J va$samp -D $diro --export=HPRC_ROOT=$root --export=HPRC_ROOT_REFERENCE=$refn --export=HPRC_ROOT_REFERENCE_HPC=$refc --export=HPRC_ROOT_REFERENCE_ODB=$odb $root/hprc/analyze.sh $samp > $diro/analysis.jid"); }
}

1;
