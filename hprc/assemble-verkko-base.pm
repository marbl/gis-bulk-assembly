
###############################################################################
#
#  This file is part of a pipeline to (help) run Verkko assemblies on
#  numerous HPRC samples.
#
#  This is a 'United States Government Work', and is released in the public
#  domain.
#
##
sub createVerkkoBase ($$$$$$$) {
  my $samp  = shift @_;
  my $flav  = "verkko-base";
  my $hifi  = shift @_;
  my $nano  = shift @_;
  my $missi = shift @_;
  my $unava = shift @_;
  my $compl = shift @_;
  my $params = shift @_;
  my $sdir  = "$rasm/$samp";

  my $ogfa = "5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip";

  if (!$missi && !$compl && !$unava && !-e "$sdir/$flav.sh") {
    system("mkdir -p $sdir")  if (! -d $sdir);

    open(CMD, "> $sdir/$flav.sh") or die "Failed to open '$sdir/$flav.sh' for writing: $!\n";
    print CMD "#!/bin/sh\n";
    print CMD "#\n";
    print CMD "#SBATCH --cpus-per-task=2\n";
    print CMD "#SBATCH --mem=16g\n";
    print CMD "#SBATCH --time=4-0\n";
    print CMD "#SBATCH --output=$sdir/$flav.%j.log\n";
    print CMD "#SBATCH --job-name=base$samp\n";
    print CMD "#\n";
    print CMD "set -o pipefail\n";
    print CMD "set -x\n";
    print CMD "\n";
    print CMD "mkdir -p $sdir\n";
    print CMD "cd       $sdir\n";
    print CMD "mkdir -p $flav\n";
    print CMD "\n";
    print CMD "export REF_CACHE=$ENV{'REF_CACHE'}\n";
    print CMD "\n";
    print CMD "module load python\n";
    print CMD "module load winnowmap\n";
    print CMD "module load mashmap\n";
    print CMD "module load samtools\n";
    print CMD "module load bwa\n";
    print CMD "module load seqtk\n";
    print CMD "\n";
    print CMD "#\n";
    print CMD "#  Run an initial generic assembly to make corrected reads and the initial graph.\n";
    print CMD "#  Also makes a contig fasta for later analysis.\n";
    print CMD "#\n";
    print CMD "\n";
    print CMD "if [ ! -e $flav/contigs.fasta ] ; then\n";
    print CMD "  if [ ! -e '$flav/emptyfile' ] ; then touch $flav/emptyfile ; fi\n";
    print CMD "  $root/$rsoft/verkko/bin/verkko --graphaligner $root/$rsoft/graphaligner/bin/GraphAligner --slurm -d $flav \\\n";
    print CMD "    --snakeopts '--until untip' \\\n";
    print CMD "    --ovb-run 8 32 32 \\\n";
    print CMD "    " . getScreenOption("$root/$rsoft/verkko") . " \\\n";
    print CMD "    --hifi $hifi \\\n";
    print CMD "    --nano $nano \\\n";
    print CMD "    $params      \\\n";
    print CMD "  && \\\n";
    print CMD "  ln -s $ogfa.gfa $flav/contigs.gfa \\\n";
    print CMD "  && \\\n";
    print CMD "  awk '/^S/ { print \">\"\$2; print \$3 }' \\\n";
    print CMD "    < $flav/contigs.gfa \\\n";
    print CMD "    | fold -c \\\n";
    print CMD "    > $flav/contigs.fasta \\\n";
    print CMD "  && \\\n";
    print CMD "  rm -rf $flav/0-correction \\\n";
    print CMD "  && \\\n";
    print CMD "  tar -C $flav -cf - .snakemake    | gzip -9 > $flav/snakemake.tar.gz     && rm -rf $flav/.snakemake && \\\n";
    print CMD "  tar -C $flav -cf - batch-scripts | gzip -9 > $flav/batch-scripts.tar.gz && rm -rf $flav/batch-scripts\n";
    print CMD "fi\n";
    print CMD "\n";
    print CMD "rm -f $sdir/$flav.jid\n";
  }

  return "$sdir/$flav";
}

1;
