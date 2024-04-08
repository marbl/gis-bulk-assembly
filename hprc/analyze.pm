
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
  my $submit = exists $$opts{"submit"};
  my $diro   = "$root/assemblies/$samp/chr-assign";

  #  Check that inputs exist.

  if (! -e "$root/assemblies/$samp/assembly.fasta") {
    print STDERR "$samp - Assembly not present.\n";
    return
  }

  #  Make a place to work.

  if (! -d "$diro") {
    system("mkdir -p $diro");
  }

  #  Run, if not finished already.

  if (-e "$diro/chr-assign.out") {
    print STDERR "$samp - Chromosome Assignment Finished (remove chr-assign/ to rerun).\n";
  }
  else {
    print STDERR "$samp - Chromosome Assignment Running.\n";

    open(CMD, "> $diro/chr-assign.sh") or die "Failed to open 'assemblies/$samp/chr-assign.sh' for writing: $!\n";

    print CMD "#!/bin/sh\n";
    print CMD "#\n";
    print CMD "#SBATCH --cpus-per-task=4\n";
    print CMD "#SBATCH --mem=4g\n";
    print CMD "#SBATCH --time=1:00:00\n";       #  Actually like 5 minutes
    print CMD "#SBATCH --output=$root/assemblies/$samp.%j.log\n";
    print CMD "#SBATCH --job-name=va$samp\n";
    print CMD "#\n";
    print CMD "set -e\n";
    print CMD "set -x\n";
    print CMD "\n";
    print CMD "mkdir -p $diro\n";
    print CMD "cd       $diro\n";
    print CMD "\n";
    print CMD "export REF_CACHE=$root/samtools-ref-cache\n";
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
    print CMD "perl $root/hprc/chromosome-assign.pl < asm-chr.mashmap.out\n";
    print CMD "\n";
    print CMD "exit 0\n";

    system("cd $diro/ && sh ./chr-assign.sh > chr-assign.out 2> chr-assign.err");

    print STDERR "$samp - Chromosome Assignment Finished.\n";
  }
}


sub startTelomereAnalysis ($$) {
  my $samp   = shift @_;
  my $opts   = shift @_;
  my $submit = exists $$opts{"submit"};
  my $diro   = "$root/assemblies/$samp/analysis";

  #  Check that inputs exist.

  if (! -e "$root/assemblies/$samp/assembly.fasta") {
    print STDERR "$samp - Assembly not present.\n";
    return
  }

  #  Make a place to work.

  if (! -d "$diro") {
    system("mkdir -p $diro");
  }

  #  Run, if not finished already.

  if (-e "$root/assemblies/$samp/analysis.out") {
    print STDERR "$samp - Telomere Analysis finished (remove analysis.out to rerun).\n";
  }
  else {
    print STDERR "$samp - Telomere Analysis Running.\n";
    system("cd $root/assemblies/$samp && sh ../../hprc/analyze.sh $samp > analysis.out");  # 2> analysis.err
    print STDERR "$samp - Telomere Analysis Finished.\n";
  }
}

1;
