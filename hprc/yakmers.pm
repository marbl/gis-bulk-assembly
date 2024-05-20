
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
use hprc::aws;

#y $yak = from module
my $seqr = "$root/software/seqrequester/build/bin/seqrequester";

sub emitYak ($$$$) {
  my $yako = shift @_;
  my $samp = shift @_;
  my $type = shift @_;
  my $inps = shift @_;

  my $p1 = "$yako/yak-$samp-$type-1.fifo";
  my $p2 = "$yako/yak-$samp-$type-2.fifo";

  if ($inps eq "") {
    print CMD "#\n";
    print CMD "#  Type $type has no data in the bucket, so nothing to compute!\n";
    print CMD "#\n";
  }
  else {
    print CMD "#\n";
    print CMD "#  Type $type.\n";
    print CMD "#\n";
  }

  print CMD "\n";
  print CMD "if [ \"$inps\" != \"\" -a \\\n";
  print CMD "     ! -e $yako/$type.yak ] ; then\n";
  print CMD "  mkfifo $p1\n";
  print CMD "  mkfifo $p2\n";
  print CMD "\n";
  print CMD "  $seqr extract -fasta $inps > $p1 &\n";
  print CMD "  $seqr extract -fasta $inps > $p2 &\n";
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
  print CMD "\n";
  print CMD "\n";
}


sub computeYakmers ($$) {
  my $samp = shift @_;
  my $opts = shift @_;
  my $subm = exists($$opts{"submit"});
  my $yako = "$root/hprc-data/$samp/yakmers";   #  Yakko waves to Wakko and Dot.

  #  Decide which yak we need to compute.  Only child, only parent, or both?

  my $cdata = (numFiles($samp, "ilmn") > 0);
  my $pdata = (numFiles($samp, "mat-ilmn") > 0) && (numFiles($samp, "pat-ilmn") > 0);

  #  Check if outputs exist.

  my $ccompl   = ($cdata == 0) || ((-e "$yako/ilmn.yak"));
  my $pcompl   = ($pdata == 0) || ((-e "$yako/mati.yak") && (-e "$yako/pati.yak"));

  my $compl    = $ccompl && $pcompl;

  #  Check that inputs exist.

  my ($nilmn, $ilmn) = (numFiles($samp, "ilmn"),          getDownloadedFiles($samp, "ilmn"));
  my ($nmati, $mati) = (numFiles($samp, "mat-ilmn"),      getDownloadedFiles($samp, "mat-ilmn"));
  my ($npati, $pati) = (numFiles($samp, "pat-ilmn"),      getDownloadedFiles($samp, "pat-ilmn"));

  my $cdownl = ((($nilmn == 0) || ($ilmn ne "")));    #  True if there are no reads in the data set or                                                                                                                     
  my $pdownl = ((($nmati == 0) || ($mati ne "")) &&   #  all reads are downloaded.  getDownloadedFiles()                                                                                                                        
                (($npati == 0) || ($pati ne "")));    #  returns empty-string if any file is missing.                                                                                                                           

  my $downl = $cdownl && $pdownl;

  #  Be a little more informative about status.

  my $cstat;
  my $pstat;

  if    ($cdata  == 0)    { $cstat = "no-data";      }   #  No data exists in bucket.
  elsif ($ccompl == 1)    { $cstat = "finished";     }   #  Computed OR no data exists.
  elsif ($cdownl == 0)    { $cstat = "not-fetched";  }   #  Not computed, data exists but not downloaded
  else                    { $cstat = "ready-to-run"; }   #  Need to compute.

  if    ($pdata  == 0)    { $pstat = "no-data";      }
  elsif ($pcompl == 1)    { $pstat = "finished";     }
  elsif ($pdownl == 0)    { $pstat = "not-fetched";  }
  else                    { $pstat = "ready-to-run"; }

  #  Launch.

  if ($downl && !$compl) {
    system("mkdir -p $yako");
    open(CMD, "> $yako/yakmers.sh") or die "Failed to open '$yako/yakmers.sh' for writing: $!\n";
    print CMD "#!/bin/sh\n";
    print CMD "#\n";
    #rint CMD "#SBATCH --gres=lscratch:512\n";
    print CMD "#SBATCH --cpus-per-task=16\n";
    print CMD "#SBATCH --mem=96g\n";
    print CMD "#SBATCH --time=8:00:00\n";
    print CMD "#SBATCH --output=$yako/yakmers.err\n";
    print CMD "#SBATCH --job-name=yak$samp\n";
    print CMD "#\n";
    print CMD "set -x\n";
    print CMD "\n";
    print CMD "export REF_CACHE=$ENV{'REF_CACHE'}\n";
    print CMD "\n";
    print CMD "module load yak\n";
    print CMD "yak version\n";
    print CMD "\n";
    print CMD "cd $yako\n";
    print CMD "\n";

    print CMD "yakCPUs=\$SLURM_JOB_CPUS_PER_NODE\n";
    print CMD "\n";
    print CMD "if [ x\$yakCPUs = x ] ; then\n";
    print CMD "  echo 'WARNING: SLURM_JOB_CPUS_PER_NODE not set, using minimal CPUs instead.'\n";
    print CMD "  yakCPUs=2\n";
    print CMD "fi\n";
    print CMD "\n";

    emitYak($yako, $samp, "ilmn", $ilmn);
    emitYak($yako, $samp, "mati", $mati);
    emitYak($yako, $samp, "pati", $pati);

    print CMD "\n";
    print CMD "\n";
    print CMD "unlink $yako/yakmers.jid\n";
    print CMD "\n";
    print CMD "exit 0\n";
    close(CMD);
  }

  if    ($compl)                      { printf "$samp/yakmers - FINISHED         (child:%-12s  parent:%s)\n", $cstat, $pstat; }
  elsif (-e "$yako/yakmers.jid")      { printf "$samp/yakmers - RUNNING          (child:%-12s  parent:%s)\n", $cstat, $pstat; }
  elsif (-e "$yako/yakmers.err")      { printf "$samp/yakmers - CRASHED          (child:%-12s  parent:%s)\n", $cstat, $pstat; }
  elsif (! $downl)                    { printf "$samp/yakmers - NOT-FETCHED      (child:%-12s  parent:%s)\n", $cstat, $pstat; }
  elsif (! $$opts{"submit"})          { printf "$samp/yakmers - READY-TO-COMPUTE (child:%-12s  parent:%s)\n", $cstat, $pstat; }
  else                                { printf "$samp/yakmers - SUBMITTED        (child:%-12s  parent:%s)\n", $cstat, $pstat; system("sbatch $yako/yakmers.sh > $yako/yakmers.jid"); }
}

1;
