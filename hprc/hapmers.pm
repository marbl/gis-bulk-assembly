
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
  my $samp = shift @_;
  my $opts = shift @_;
  my $hapo   = "$root/aws-data/$samp/hapmers";

  print "\n";
  print "$samp\n";

  #  Check if outputs exist.

  if ((-e "$hapo/mati.hapmers.meryl/merylIndex") &&
      (-e "$hapo/pati.hapmers.meryl/merylIndex")) {
    print "  Outputs $samp/mati.hapmers.meryl/ and\n";
    print "          $samp/pati.hapmers.meryl/ exist, nothing to do.\n";
    return;
  }

  #  Make a place to work.

  if (! -d "$hapo") {
    system("mkdir -p $hapo");
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

  if (-e "$hapo/hapmers.jid") {
    print "  Currently running.\n";
    return;
  }

  my $meryl = "$root/software/meryl/build/bin/meryl";
  my $d     = "/lscratch/\$SLURM_JOB_ID";

  #  Launch.
  #    HG04199 used just over 250GB disk.

  open(CMD, "> $hapo/hapmers.sh") or die "Failed to open '$hapo/hapmers.sh' for writing: $!\n";
  print CMD "#!/bin/sh\n";
  print CMD "#\n";
  print CMD "#SBATCH --gres=lscratch:512\n";
  print CMD "#SBATCH --cpus-per-task=16\n";
  print CMD "#SBATCH --mem=128g\n";
  print CMD "#SBATCH --time=1-0\n";
  print CMD "#SBATCH --output=$hapo/hapmers.err\n";
  print CMD "#SBATCH --job-name=hapmer$samp\n";
  print CMD "#\n";
  print CMD "set -x\n";
  print CMD "\n";
  print CMD "export REF_CACHE=$root/samtools-ref-cache\n";
  print CMD "\n";
  print CMD "cd $hapo\n";
  print CMD "\n";
  print CMD "\n";
  print CMD "if [ ! -e $d/ilmn.meryl/merylIndex ] ; then\n";
  print CMD "  $meryl k=31 threads=16 memory=96g count compress \\\n";
  print CMD "    output $d/ilmn.meryl \\\n";
  print CMD "         $ilmn \\\n";
  print CMD "  || \\\n";
  print CMD "  rm -rf $d/ilmn.meryl\n";
  print CMD "fi\n";
  print CMD "\n";
  print CMD "if [ ! -e $d/mati.meryl/merylIndex ] ; then\n";
  print CMD "  $meryl k=31 threads=16 memory=96g count compress \\\n";
  print CMD "    output $d/mati.meryl \\\n";
  print CMD "         $mati \\\n";
  print CMD "  || \\\n";
  print CMD "  rm -rf $d/mati.meryl\n";
  print CMD "fi\n";
  print CMD "\n";
  print CMD "if [ ! -e $d/pati.meryl/merylIndex ] ; then\n";
  print CMD "  $meryl k=31 threads=16 memory=96g count compress \\\n";
  print CMD "    output $d/pati.meryl \\\n";
  print CMD "         $pati \\\n";
  print CMD "  || \\\n";
  print CMD "  rm -rf $d/pati.meryl\n";
  print CMD "fi\n";
  print CMD "\n";
  print CMD "\n";

  my $filt1 = "\$(java -jar -Xmx1g $root/software/merqury/eval/kmerHistToPloidyDepth.jar";
  my $filt2 = "2> /dev/null | awk '\$1 == 0 { print \$3 }')";

  print CMD "#\n";
  print CMD "#  Automagically generated.\n";
  print CMD "#  Vastly more readable in hprc/hapmers.pm.\n";
  print CMD "#\n";
  print CMD "\n";
  print CMD "#  maternal/paternal specific kmers\n";
  print CMD "$meryl difference $d/mati.meryl $d/pati.meryl output $d/mati.only.meryl || rm -rf $d/mati.only.meryl\n";
  print CMD "$meryl difference $d/pati.meryl $d/mati.meryl output $d/pati.only.meryl || rm -rf $d/pati.only.meryl\n";
  print CMD "\n";
  print CMD "$meryl histogram $d/mati.only.meryl > $d/mati.only.histogram || rm -f $d/mati.only.histogram\n";
  print CMD "$meryl histogram $d/pati.only.meryl > $d/pati.only.histogram || rm -f $d/pati.only.histogram\n";
  print CMD "\n";
  print CMD "tm=$filt1 $d/mati.only.histogram $filt2\n";
  print CMD "tp=$filt1 $d/pati.only.histogram $filt2\n";
  print CMD "\n";
  print CMD "echo Filter hap-specific mat/pat gt \$tm/\$tp\n";
  print CMD "\n";
  print CMD "$meryl intersect $d/ilmn.meryl [ greater-than \$tm $d/mati.only.meryl ] output $d/mati.inherited.meryl || rm -rf $d/mati.inherited.meryl\n";
  print CMD "$meryl intersect $d/ilmn.meryl [ greater-than \$tp $d/pati.only.meryl ] output $d/pati.inherited.meryl || rm -rf $d/pati.inherited.meryl\n";
  print CMD "\n";
  print CMD "$meryl histogram $d/mati.inherited.meryl > $d/mati.inherited.histogram || rm -f $d/mati.inherited.histogram\n";
  print CMD "$meryl histogram $d/pati.inherited.meryl > $d/pati.inherited.histogram || rm -f $d/pati.inherited.histogram\n";
  print CMD "\n";
  print CMD "tm=$filt1 $d/mati.inherited.histogram $filt2\n";
  print CMD "tp=$filt1 $d/pati.inherited.histogram $filt2\n";
  print CMD "\n";
  print CMD "echo Filter inherited mat/pat gt \$tm/\$tp\n";
  print CMD "\n";
  print CMD "$meryl greater-than \$tm $d/mati.inherited.meryl output $d/mati.hapmers.meryl || rm -rf $d/mati.hapmers.meryl\n";
  print CMD "$meryl greater-than \$tp $d/pati.inherited.meryl output $d/pati.hapmers.meryl || rm -rf $d/pati.hapmers.meryl\n";
  print CMD "\n";
  print CMD "\n";
  print CMD "#sync -av $d/ilmn.meryl .\n";            #  About 50GB each
  print CMD "#sync -av $d/mati.meryl .\n";
  print CMD "#sync -av $d/pati.meryl .\n";
  print CMD "rsync -av $d/mati.hapmers.meryl .\n";   #  2 GB each.
  print CMD "rsync -av $d/pati.hapmers.meryl .\n";
  print CMD "\n";
  print CMD "rm $hapo/hapmers.jid\n";
  print CMD "\n";
  print CMD "exit 0\n";
  close(CMD);

  print STDOUT "sbatch $hapo/hapmers.sh > $hapo/hapmers.jid\n";
  system("sbatch $hapo/hapmers.sh > $hapo/hapmers.jid 2>&1")   if (exists($$opts{"submit"}));
}

1;
