
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
use hprc::aws;

my $meryl = "$root/software/meryl/build/bin/meryl";

sub computeHapmers ($$) {
  my $samp = shift @_;
  my $opts = shift @_;
  my $hapo   = "$root/hprc-data/$samp/hapmers";
  my $locd  = "/lscratch/\$SLURM_JOB_ID";

  #  Check if outputs exist.

  my $compl = ((-e "$hapo/mati.hapmers.meryl/merylIndex") &&
               (-e "$hapo/pati.hapmers.meryl/merylIndex"));

  #  Check that inputs exist.

  my ($nilmn, $ilmn) = (numFiles($samp, "ilmn"),          getDownloadedFiles($samp, "ilmn"));
  my ($nmati, $mati) = (numFiles($samp, "mat-ilmn"),      getDownloadedFiles($samp, "mat-ilmn"));
  my ($npati, $pati) = (numFiles($samp, "pat-ilmn"),      getDownloadedFiles($samp, "pat-ilmn"));

  my $downl = ((($nilmn == 0) || ($ilmn ne "")) &&   #  True if there are no reads in the data set or                                                                                                                     
               (($nmati == 0) || ($mati ne "")) &&   #  all reads are downloaded.  getDownloadedFiles()                                                                                                                        
               (($npati == 0) || ($pati ne "")));    #  returns empty-string if any file is missing.                                                                                                                           

  #  Create a script - if data is donwloaded and script doesn't exist.                                                                                                                                                         

  if ($downl && !$compl) {
    system("mkdir -p $hapo")   if (! -d "$hapo");
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
    print CMD "export REF_CACHE=$ENV{'REF_CACHE'}\n";
    print CMD "\n";
    print CMD "cd $hapo\n";
    print CMD "\n";
    print CMD "\n";
    print CMD "if [ ! -e $locd/ilmn.meryl/merylIndex ] ; then\n";
    print CMD "  $meryl k=31 threads=16 memory=96g count compress \\\n";
    print CMD "    output $locd/ilmn.meryl \\\n";
    print CMD "         $ilmn \\\n";
    print CMD "  || \\\n";
    print CMD "  rm -rf $locd/ilmn.meryl\n";
    print CMD "fi\n";
    print CMD "\n";
    print CMD "if [ ! -e $locd/mati.meryl/merylIndex ] ; then\n";
    print CMD "  $meryl k=31 threads=16 memory=96g count compress \\\n";
    print CMD "    output $locd/mati.meryl \\\n";
    print CMD "         $mati \\\n";
    print CMD "  || \\\n";
    print CMD "  rm -rf $locd/mati.meryl\n";
    print CMD "fi\n";
    print CMD "\n";
    print CMD "if [ ! -e $locd/pati.meryl/merylIndex ] ; then\n";
    print CMD "  $meryl k=31 threads=16 memory=96g count compress \\\n";
    print CMD "    output $locd/pati.meryl \\\n";
    print CMD "         $pati \\\n";
    print CMD "  || \\\n";
    print CMD "  rm -rf $locd/pati.meryl\n";
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
    print CMD "$meryl difference $locd/mati.meryl $locd/pati.meryl output $locd/mati.only.meryl || rm -rf $locd/mati.only.meryl\n";
    print CMD "$meryl difference $locd/pati.meryl $locd/mati.meryl output $locd/pati.only.meryl || rm -rf $locd/pati.only.meryl\n";
    print CMD "\n";
    print CMD "$meryl histogram $locd/mati.only.meryl > $locd/mati.only.histogram || rm -f $locd/mati.only.histogram\n";
    print CMD "$meryl histogram $locd/pati.only.meryl > $locd/pati.only.histogram || rm -f $locd/pati.only.histogram\n";
    print CMD "\n";
    print CMD "tm=$filt1 $locd/mati.only.histogram $filt2\n";
    print CMD "tp=$filt1 $locd/pati.only.histogram $filt2\n";
    print CMD "\n";
    print CMD "echo Filter hap-specific mat/pat gt \$tm/\$tp\n";
    print CMD "\n";
    print CMD "$meryl intersect $locd/ilmn.meryl [ greater-than \$tm $locd/mati.only.meryl ] output $locd/mati.inherited.meryl || rm -rf $locd/mati.inherited.meryl\n";
    print CMD "$meryl intersect $locd/ilmn.meryl [ greater-than \$tp $locd/pati.only.meryl ] output $locd/pati.inherited.meryl || rm -rf $locd/pati.inherited.meryl\n";
    print CMD "\n";
    print CMD "$meryl histogram $locd/mati.inherited.meryl > $locd/mati.inherited.histogram || rm -f $locd/mati.inherited.histogram\n";
    print CMD "$meryl histogram $locd/pati.inherited.meryl > $locd/pati.inherited.histogram || rm -f $locd/pati.inherited.histogram\n";
    print CMD "\n";
    print CMD "tm=$filt1 $locd/mati.inherited.histogram $filt2\n";
    print CMD "tp=$filt1 $locd/pati.inherited.histogram $filt2\n";
    print CMD "\n";
    print CMD "echo Filter inherited mat/pat gt \$tm/\$tp\n";
    print CMD "\n";
    print CMD "$meryl greater-than \$tm $locd/mati.inherited.meryl output $locd/mati.hapmers.meryl || rm -rf $locd/mati.hapmers.meryl\n";
    print CMD "$meryl greater-than \$tp $locd/pati.inherited.meryl output $locd/pati.hapmers.meryl || rm -rf $locd/pati.hapmers.meryl\n";
    print CMD "\n";
    print CMD "\n";
    print CMD "#sync -av $locd/ilmn.meryl .\n";           #  About 50GB each; we don't need them later.
    print CMD "#sync -av $locd/mati.meryl .\n";
    print CMD "#sync -av $locd/pati.meryl .\n";
    print CMD "rsync -av $locd/mati.hapmers.meryl .\n";   #  2 GB each.
    print CMD "rsync -av $locd/pati.hapmers.meryl .\n";
    print CMD "\n";
    print CMD "rm $hapo/hapmers.jid\n";
    print CMD "\n";
    print CMD "exit 0\n";
    close(CMD);
  }

  if    ($compl)                      { print "$samp/hapmers - FINISHED\n"; }
  elsif (-e "$hapo/hapmers.jid")      { print "$samp/hapmers - RUNNING\n"; }
  elsif (-e "$hapo/hapmers.err")      { print "$samp/hapmers - CRASHED\n"; }
  elsif (! $downl)                    { print "$samp/hapmers - NOT-FETCHED\n"; }
  elsif (! $$opts{"submit"})          { print "$samp/hapmers - READY-TO-COMPUTE\n"; }
  else                                { print "$samp/hapmers - SUBMITTED\n"; system("sbatch $hapo/hapmers.sh > $hapo/hapmers.jid"); }
}

1;
