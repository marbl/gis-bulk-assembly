
###############################################################################
#
#  This file is part of a pipeline to (help) run Verkko assemblies on
#  numerous HPRC samples.
#
#  This is a 'United States Government Work', and is released in the public
#  domain.
#
##


sub locateHapmers ($) {     #  In array context, returns the expected locations of the two hapmer
  my $samp  = shift @_;     #  databases.  In scalar context, returns false/true if those exist.
  my $mati  = "$root/hprc-data/$samp/hapmers/mati.hapmers.meryl";
  my $pati  = "$root/hprc-data/$samp/hapmers/pati.hapmers.meryl";
  my $missi = ((! -e "$mati/merylIndex") ||
               (! -e "$pati/merylIndex"));

  return wantarray() ? ($mati, $pati) : $missi;
}

sub createVerkkoTrio ($$$$$$) {
  my $samp  = shift @_;
  my $flav  = "verkko-trio";
  my $hifi  = shift @_;
  my $nano  = shift @_;
  my $missi = shift @_;
  my $unava = shift @_;
  my $compl = shift @_;
  my $sdir  = "$rasm/$samp";

  my ($mati, $pati) = locateHapmers($samp);

  if (!$missi && !$compl && !$unava) {
    system("mkdir -p $sdir")  if (! -d $sdir);

    open(CMD, "> $sdir/$flav.sh") or die "Failed to open '$sdir/$flav.sh' for writing: $!\n";
    print CMD "#!/bin/sh\n";
    print CMD "#\n";
    print CMD "#SBATCH --cpus-per-task=2\n";
    print CMD "#SBATCH --mem=16g\n";
    print CMD "#SBATCH --time=4-0\n";
    print CMD "#SBATCH --output=$sdir/$flav.%j.log\n";
    print CMD "#SBATCH --job-name=trio$samp\n";
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
    print CMD "module load winnowmap\n";
    print CMD "module load mashmap\n";
    print CMD "module load samtools\n";
    print CMD "module load bwa\n";
    print CMD "module load seqtk\n";
    print CMD "\n";
    print CMD "#\n";
    print CMD "#  Fail if no base assembly.\n";
    print CMD "#\n";
    print CMD "if [ ! -e verkko-base/contigs.fasta ] ; then\n";
    print CMD "  echo \"No base assembly found, can't start trio assembly.\"\n";
    print CMD "\n";
    print CMD "#\n";
    print CMD "#  Continue a base assembly into a trio assembly.\n";
    print CMD "#\n";
    print CMD "elif [ ! -e $flav/assembly.fasta ] ; then\n";
    print CMD "  if [ ! -e '$flav/emptyfile' ] ;               then  touch -r verkko-base/emptyfile               $flav/emptyfile ;               fi\n";
    print CMD "  if [ ! -e '$flav/1-buildGraph' ] ;            then  ln -s ../verkko-base/1-buildGraph            $flav/1-buildGraph ;            fi\n";
    print CMD "  if [ ! -e '$flav/2-processGraph' ] ;          then  ln -s ../verkko-base/2-processGraph          $flav/2-processGraph ;          fi\n";
    print CMD "  if [ ! -e '$flav/3-align' ] ;                 then  ln -s ../verkko-base/3-align                 $flav/3-align ;                 fi\n";
    print CMD "  if [ ! -e '$flav/3-alignTips' ] ;             then  ln -s ../verkko-base/3-alignTips             $flav/3-alignTips ;             fi\n";
    print CMD "  if [ ! -e '$flav/4-processONT' ] ;            then  ln -s ../verkko-base/4-processONT            $flav/4-processONT ;            fi\n";
    print CMD "  if [ ! -e '$flav/5-untip' ] ;                 then  ln -s ../verkko-base/5-untip                 $flav/5-untip ;                 fi\n";
    print CMD "  if [ ! -e '$flav/hifi-corrected.fasta.gz' ] ; then  ln -s ../verkko-base/hifi-corrected.fasta.gz $flav/hifi-corrected.fasta.gz ; fi\n";
    print CMD "\n";
    print CMD "  $root/software/verkko/bin/verkko --graphaligner $root/software/graphaligner/bin/GraphAligner --slurm -d $flav \\\n";
    print CMD "    --ovb-run 8 32 32 \\\n";
    print CMD "    --screen human \\\n";
    print CMD "    --hifi $hifi \\\n";
    print CMD "    --nano $nano \\\n";
    print CMD "    --hap-kmers $mati \\\n";
    print CMD "                $pati \\\n";
    print CMD "                trio \\\n";
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
