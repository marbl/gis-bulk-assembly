
###############################################################################
#
#  This file is part of a pipeline to (help) run Verkko assemblies on
#  numerous HPRC samples.
#
#  This is a 'United States Government Work', and is released in the public
#  domain.
#
##

#
#  Archive intermediate files.  Assumes that cleanupAssembly() has been
#  performed, but you can happily waste disk and archive the full outputs.
#

sub archiveAssembly($$) {
  my $samp = shift @_;
  my $opts = shift @_;
  my @flavs;
  #y $flav = $$opts{"flavor"};
  my $subm = $$opts{"submit"};

  if ($$opts{"flavor"} ne "") {
    push @flavs, $$opts{"flavor"};
  }
  else {
    push @flavs, "verkko-trio";
    push @flavs, "verkko-hi-c";
    push @flavs, "verkko-full";
    push @flavs, "verkko";
  }

  foreach my $flav (@flavs) {
    my $path = "$rasm/$samp/$flav";

    next   if (!isFinished($samp, $flav));

    foreach my $adir (qw(1-buildGraph 2-processGraph 3-align 3-alignTips 4-processONT 5-untip 6-layoutContigs 7-consensus 8-hicPipeline)) {
      if ((-e "$path/$adir") && (! -e "$path/$adir.tar.gz")) {
        #y $samp = shift @_;   #  Sample name, for logging.
        #y $flav = shift @_;   #  Flavor, for logging.
        #y $path = shift @_;   #  Path to directory containing directory to archive.
        #y $adir = shift @_;   #  Name of directory to archive.

        open(CMD, "> ./$adir.tar.sh") or die "Failed to open '$adir.tar.sh' for writing: $!\n";
        print CMD "#!/bin/sh\n";
        print CMD "#\n";
        print CMD "#SBATCH --cpus-per-task=12\n";
        print CMD "#SBATCH --mem=16g\n";
        print CMD "#SBATCH --time=4-0\n";
        print CMD "#SBATCH --output=./$adir.tar.gz.err\n";
        print CMD "#SBATCH --job-name=gz$samp$flav$adir\n";
        print CMD "#\n";
        print CMD "set -o pipefail\n";
        print CMD "set -x\n";
        print CMD "\n";
        print CMD "module load pigz\n";
        print CMD "\n";
        print CMD "cd $path/\n";
        print CMD "tar -cf - ./$adir | pigz -9 -p 12 -c > ./$adir.tar.gz \\\n";
        print CMD "\n";
        print CMD "if [ \$? = 0 ] ; then\n";
        print CMD "  chmod 444 ./$adir.tar.gz\n";
        print CMD "  chmod -R +w ./$adir\n";
        print CMD "  rm -rf ./$adir ./$adir.tar.sh ./$adir.tar.jid ./$adir.tar.gz.err\n";
        print CMD "  exit 0\n";
        print CMD "else\n";
        print CMD "  rm -f ./$adir.tar.gz ./$adir.tar.jid\n";
        print CMD "  exit 1\n";
        print CMD "fi\n";
        close(CMD);

        print STDOUT "  sbatch $path/$adir.tar.sh > $path/$adir.tar.jid", (($subm) ? "" : " (not submitted)"), "\n";
        system("sbatch $path/$adir.tar.sh > $path/$adir.tar.jid 2>&1")   if ($subm);
      }
    }
  }
}

1;
