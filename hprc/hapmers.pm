
###############################################################################
#
#  This file is part of a pipeline to (help) run Verkko assemblies on
#  numerous HPRC samples.
#
#  This is a 'United States Government Work', and is released in the public
#  domain.
#
##

package hprc::hapmers;
require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(computeHapmers);

use strict;
use warnings "all";
no  warnings "uninitialized";

use hprc::samples;
use hprc::assemble;


sub computeHapmers ($$) {
  my $samp   = shift @_;
  my $submit = shift @_;

  #  Check if outputs exist.

  if ((-e "$root/assemblies/$samp/mati.inherited.meryl/merylIndex") &&
      (-e "$root/assemblies/$samp/pati.inherited.meryl/merylIndex")) {
    print "Outputs $samp/mati.inherited.meryl/ and\n";
    print "        $samp/pati.inherited.meryl/ exist, nothing to do.\n";
    return;
  }

  #  Make a place to work.

  if (! -d "$root/assemblies/$samp") {
    system("mkdir -p $root/assemblies/$samp");
  }

  #  Check that inputs exist.

  my $ilmn = getFiles($samp, "ilmn");
  my $mati = getFiles($samp, "mat-ilmn");
  my $pati = getFiles($samp, "pat-ilmn");

  #  Make sure it isn't running.

  #  Launch.
  #    HG04199 used just over 250GB disk.

  open(CMD, "> $root/assemblies/$samp-hapmers.sh") or die "Failed to open 'assemblies/$samp-hapmers.sh' for writing: $!\n";
  print CMD "#!/bin/sh\n";
  print CMD "#\n";
  print CMD "#SBATCH --gres=lscratch:512\n";
  print CMD "#SBATCH --cpus-per-task=16\n";
  print CMD "#SBATCH --mem=128g\n";
  print CMD "#SBATCH --time=1-0\n";
  print CMD "#SBATCH --output=$root/assemblies/$samp-hapmers.%j.log\n";
  print CMD "#SBATCH --job-name=hapmer$samp\n";
  print CMD "#\n";
  print CMD "set -e\n";
  print CMD "set -x\n";
  print CMD "\n";
  print CMD "export REF_CACHE=$root/samtools-ref-cache\n";
  print CMD "\n";
  print CMD "cd $root/assemblies/$samp\n";
  print CMD "\n";
  print CMD "\n";
  print CMD "if [ ! -e /lscratch/\$SLURM_JOB_ID/ilmn.meryl/merylIndex ] ; then\n";
  print CMD "  $root/software/meryl/build/bin/meryl \\\n";
  print CMD "    k=31 threads=16 memory=96g \\\n";
  print CMD "    count \\\n";
  print CMD "      output /lscratch/\$SLURM_JOB_ID/ilmn.meryl \\\n";
  print CMD "        compress \\\n";
  print CMD "         $ilmn\n";
  print CMD "fi\n";
  print CMD "\n";
  print CMD "\n";
  print CMD "\n";
  print CMD "if [ ! -e /lscratch/\$SLURM_JOB_ID/mati.meryl/merylIndex ] ; then\n";
  print CMD "  $root/software/meryl/build/bin/meryl \\\n";
  print CMD "    k=31 threads=16 memory=96g \\\n";
  print CMD "    count \\\n";
  print CMD "      output /lscratch/\$SLURM_JOB_ID/mati.meryl \\\n";
  print CMD "        compress \\\n";
  print CMD "         $mati\n";
  print CMD "fi\n";
  print CMD "\n";
  print CMD "\n";
  print CMD "\n";
  print CMD "if [ ! -e /lscratch/\$SLURM_JOB_ID/pati.meryl/merylIndex ] ; then\n";
  print CMD "  $root/software/meryl/build/bin/meryl \\\n";
  print CMD "    k=31 threads=16 memory=96g \\\n";
  print CMD "    count \\\n";
  print CMD "      output /lscratch/\$SLURM_JOB_ID/pati.meryl \\\n";
  print CMD "        compress \\\n";
  print CMD "         $pati\n";
  print CMD "fi\n";
  print CMD "\n";
  print CMD "\n";
  print CMD "\n";

  my $meryl = "$root/software/meryl/build/bin/meryl";
  my $filt1 = "\$(java -jar -Xmx1g $root/software/merqury/eval/kmerHistToPloidyDepth.jar";
  my $filt2 = "2> /dev/null | awk '\$1 == 0 { print \$3 }')";
  my $d     = "/lscratch/\$SLURM_JOB_ID";

  #  SLACK, #meryl, Thr 7 Mar 2024, 13:09
  #
  #    @Dmitry Antipov proposed to do the intersection first; like
  #      meryl intersect $child.meryl $mat.meryl output $mat-inherited.meryl
  #    Then doing the subtraction and filtering:
  #      meryl greater-than <tresh> [ difference $mat-inherited.meryl $pat-inherited.meryl ]
  #    I wanted to try this out and see if it gives reasonable results
  #

  print CMD "#\n";
  print CMD "#  Automagically generated.\n";
  print CMD "#  Vastly more readable in hprc/hapmers.pm.\n";
  print CMD "#\n";
  print CMD "\n";
  print CMD "#  maternal specific kmers\n";
  print CMD "$meryl difference output $d/mati.only.meryl $d/mati.meryl $d/pati.meryl\n";
  print CMD "$meryl histogram $d/mati.only.meryl > $d/mati.only.histogram\n";
  print CMD "\n";
  print CMD "t=$filt1 $d/mati.only.histogram $filt2\n";
  print CMD "$meryl greater-than \$t output $d/mati.onlyfilt.meryl $d/mati.only.meryl\n";
  print CMD "\n";
  print CMD "#  paternal specific kmers\n";
  print CMD "$meryl difference output $d/pati.only.meryl $d/pati.meryl $d/mati.meryl\n";
  print CMD "$meryl histogram $d/pati.only.meryl > $d/pati.only.histogram\n";
  print CMD "\n";
  print CMD "t=$filt1 $d/pati.only.histogram $filt2\n";
  print CMD "$meryl greater-than \$t output $d/pati.onlyfilt.meryl $d/pati.only.meryl\n";
  print CMD "\n";
  print CMD "#  mat hapmers\n";
  print CMD "$meryl intersect output $d/mati.inherited.meryl $d/ilmn.meryl $d/mati.only.meryl\n";
  print CMD "$meryl histogram $d/mati.inherited.meryl > $d/mati.inherited.histogram\n";
  print CMD "\n";
  print CMD "t=$filt1 $d/mati.inherited.histogram $filt2\n";
  print CMD "$meryl greater-than \$t output $d/mati.hapmer $d/mati.inherited.meryl\n";
  print CMD "\n";
  print CMD "#  pat hapmers\n";
  print CMD "$meryl intersect output $d/pati.inherited.meryl $d/ilmn.meryl $d/pati.only.meryl \n";
  print CMD "$meryl histogram $d/pati.inherited.meryl > $d/pati.inherited.histogram\n";
  print CMD "\n";
  print CMD "t=$filt1 $d/pati.inherited.histogram $filt2\n";
  print CMD "$meryl greater-than \$t output $d/mati.hapmer $d/mati.inherited.meryl\n";
  print CMD "\n";
  print CMD "\n";
  print CMD "\n";
  print CMD "#rsync -av $d/ilmn.meryl .\n";            #  About 50GB each
  print CMD "#rsync -av $d/mati.meryl .\n";
  print CMD "#rsync -av $d/pati.meryl .\n";
  print CMD "rsync -av $d/mati.inherited.meryl .\n";   #  2 GB each.
  print CMD "rsync -av $d/pati.inherited.meryl .\n";
  print CMD "\n";
  print CMD "exit 0\n";
  close(CMD);

  print STDOUT "sbatch $root/assemblies/$samp-hapmers.sh > $root/assemblies/$samp-hapmer.jid\n";
  system("sbatch $root/assemblies/$samp-hapmers.sh > $root/assemblies/$samp-hapmer.jid 2>&1")   if ($submit);
}

1;
