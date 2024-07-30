
###############################################################################
#
#  This file is part of a pipeline to (help) run Verkko assemblies on
#  numerous HPRC samples.
#
#  This is a 'United States Government Work', and is released in the public
#  domain.
#
##

package hprc::ribotin;
require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(startRibotinAnalysis);

use strict;
use warnings "all";
no  warnings "uninitialized";

use hprc::samples;
use hprc::assemble;


sub startRibotinAnalysis($$) {
  my $samp   = shift @_;
  my $opts   = shift @_;
  my $flav   = $$opts{"flavor"};
  my $submit = $$opts{"submit"};
  my $diro   = "$rasm/$samp/$flav/ribotin";

  #  Check if outputs exist.

  my $finished = -e "$diro/ribotin.out";

  #  Check that inputs exist.

  my $ready   = isFinished($samp, $flav);
  my $unavail = dataAvailable($samp, $flav);

  #  Create a script (we don't need to since it exists elsewhere).

  if ($ready && !$finished) {
    system("mkdir -p $diro")   if (! -d "$diro");

    open(CMD, "> $diro/ribotin.sh") or die "Failed to open '$diro/ribotin.sh' for writing: $!\n";
    print CMD "#!/bin/sh\n";
    print CMD "#\n";
    print CMD "#SBATCH --cpus-per-task=8\n";
    print CMD "#SBATCH --mem=20g\n";
    print CMD "#SBATCH --time=10:00:00\n";
    print CMD "#SBATCH --output=$diro/ribotin.%j.err\n";
    print CMD "#SBATCH --job-name=rib$samp\n";
    print CMD "#\n";
    print CMD "set -e\n";
    print CMD "set -x\n";
    print CMD "\n";
    print CMD "mkdir -p $diro\n";
    print CMD "cd       $diro\n";
    print CMD "\n";
    print CMD "module load lbzip2\n";
    print CMD "module load samtools\n";
    print CMD "module load seqtk\n";
    print CMD "\n";
    print CMD "\n";
    print CMD "if [ ! -e ../assembly.homopolymer-compressed.gfa ]; then\n";
    print CMD "   ln -s 5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.gfa ../assembly.homopolymer-compressed.gfa\n";
    print CMD "fi\n";
    print CMD "\n";
    print CMD "if [ ! -e ribotin.out ] ; then\n";
    print CMD "export PATH=\$PATH:/$root/software/ribotin/bin\n";
    # we make a copy of the ribotin template sequences because liftoff was crashing when too many DBs were pointed to the same gff
    print CMD "cp -r $root/software/ribotin/template_seqs ./\n";
    print CMD "ribotin-verkko --graphaligner $root/software/graphaligner/bin/GraphAligner --mbg $root/software/verkko/lib/verkko/bin/MBG \\\n";
    print CMD "    -t \$SLURM_CPUS_PER_TASK \\\n";
    print CMD "    -x human \\\n";
    print CMD "    --guess-tangles-using-reference ./template_seqs/chm13_rDNAs.fa \\\n";
    print CMD "    --orient-by-reference           ./template_seqs/rDNA_one_unit.fasta \\\n";
    print CMD "    --annotation-reference-fasta    ./template_seqs/rDNA_one_unit.fasta \\\n";
    print CMD "    --annotation-gff3               ./template_seqs/rDNA_annotation.gff3 \\\n";
    print CMD "    -i ../ \\\n";
    print CMD "    -o results \\\n";
    print CMD " > ribotin.WORKING \\\n";
    print CMD "&& \\\n";
    print CMD "mv ribotin.WORKING ribotin.out \\\n";
    print CMD "|| \\\n";
    print CMD "mv ribotin.WORKING ribotin.FAIL\n";
    print CMD "fi\n";
    print CMD "\n";
    print CMD "rm -f ./ribotin.jid\n";
    print CMD "rm -rf ./tmp\n";
    print CMD "rm -rf ./template_seqs\n";
    print CMD "\n";
    print CMD "exit 0\n";
  }

  #  Run, if not finished already.

  if    ($finished)                   { print "$samp/$flav/ribotin - FINISHED\n"; }
  elsif (-e "$diro/ribotin.jid")      { print "$samp/$flav/ribotin - RUNNING\n"; }
  elsif (-e "$diro/ribotin.FAIL")     { print "$samp/$flav/ribotin - CRASHED\n"; }
  elsif ($unavail)                    { print "$samp/$flav/ribotin - UNAVAILABLE ($unavail)\n"; }
  elsif (! $ready)                    { print "$samp/$flav/ribotin - ASSEMBLY-NOT-READY\n"; }
  elsif (! $$opts{"submit"})          { print "$samp/$flav/ribotin - READY-TO-COMPUTE\n"; }
  else                                { print "$samp/$flav/ribotin - SUBMITTED\n"; system("sbatch $diro/ribotin.sh > $diro/ribotin.jid"); }
}


1;
