
###############################################################################
#
#  This file is part of a pipeline to (help) run Verkko assemblies on
#  numerous HPRC samples.
#
#  This is a 'United States Government Work', and is released in the public
#  domain.
#
##


sub createVerkkoTrioHiC ($$$$$$$) {
  my $samp  = shift @_;
  my $flav  = "verkko-thic";
  my $hifi  = shift @_;
  my $nano  = shift @_;
  my $hic1  = shift @_;
  my $hic2  = shift @_;
  my $missi = shift @_;
  my $unava = shift @_;
  my $compl = shift @_;
  my $sdir  = "$rasm/$samp";

  if (!$missi && !$compl && !$unava) {
    system("mkdir -p $sdir")  if (! -d $sdir);

    open(CMD, "> $sdir/$flav.sh") or die "Failed to open '$sdir/$flav.sh' for writing: $!\n";
    print CMD "#!/bin/sh\n";
    print CMD "#\n";
    print CMD "#SBATCH --cpus-per-task=2\n";
    print CMD "#SBATCH --mem=16g\n";
    print CMD "#SBATCH --time=4-0\n";
    print CMD "#SBATCH --output=$sdir/$flav.%j.log\n";
    print CMD "#SBATCH --job-name=thic$samp\n";
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
    print CMD "#  Fail if no base or trio assembly.\n";
    print CMD "#\n";
    print CMD "if [ ! -e verkko-base/contigs.fasta ] ; then\n";
    print CMD "  echo \"No base assembly found, can't start trio+hi-c assembly.\"\n";
    print CMD "\n";
    print CMD "elif [ ! -e verkko-trio/assembly.fasta ] ; then\n";
    print CMD "  echo \"No trio assembly found, can't start trio+hi-c assembly.\"\n";
    print CMD "\n";
    print CMD "#\n";
    print CMD "#  Continue a base assembly (and a few pieces from a trio assembly) into a trio+hi-c assembly.\n";
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
    print CMD "  if [ ! -e $flav/7-consensus ] ; then\n";
    print CMD "    mkdir -p $flav/7-consensus\n";
    print CMD "    ln -s ../../verkko-trio/7-consensus/ont_subset.fasta.gz  $flav/7-consensus/ont_subset.fasta.gz\n";
    print CMD "    ln -s ../../verkko-trio/7-consensus/ont_subset.id        $flav/7-consensus/ont_subset.id\n";
    print CMD "  fi\n";
    print CMD "\n";
    print CMD "  if [ ! -e $flav/8-hicPipeline ] ; then\n";
    print CMD "    mkdir -p $flav/8-hicPipeline\n";
    print CMD "    ln -s ../../verkko-trio/6-rukki/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.noseq.gfa  $flav/8-hicPipeline/unitigs.hpc.noseq.gfa\n";
    print CMD "    ln -s ../../verkko-trio/6-rukki/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.fasta      $flav/8-hicPipeline/unitigs.hpc.fasta\n";
    print CMD "  fi\n";
    print CMD "\n";
    print CMD "  if [ ! -e $flav/8-hicPipeline/rukki.paths.tsv ] ; then\n";
    print CMD "    $root/software/verkko/bin/verkko --graphaligner $root/software/graphaligner/bin/GraphAligner --slurm -d $flav --snakeopts '-U hicPhasing' \\\n";
    print CMD "      --ovb-run 8 32 32 \\\n";
    print CMD "      --screen human \\\n";
    #rint CMD "      --rdna-scaff \\\n";   #  Verkko 2.0 needs this.
    print CMD "      --hifi $hifi \\\n";
    print CMD "      --nano $nano \\\n";
    print CMD "      --hic1 $hic1 \\\n";
    print CMD "      --hic2 $hic2 \\\n";
    print CMD "    && \\\n";
    print CMD "    mv $flav/8-hicPipeline/hicverkko.colors.tsv \\\n";
    print CMD "       $flav/8-hicPipeline/hicverkko.hiccolors.tsv \\\n";
    print CMD "    && \\\n";
    print CMD "    cp verkko-trio/6-rukki/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.colors.csv \\\n";
    print CMD "       $flav/8-hicPipeline/hicverkko.colors.tsv \\\n";
    print CMD "    && \\\n";
    print CMD "    rm -f $flav/8-hicPipeline/label* && \\\n";
    print CMD "    cp -p verkko-trio/6-rukki/label* $flav/8-hicPipeline/ \\\n";
    print CMD "    && \\\n";
    print CMD "    rm -f $flav/8-hicPipeline/rukki.paths.gaf && \\\n";
    print CMD "    rm -f $flav/8-hicPipeline/rukki.paths.tsv && \\\n";
    print CMD "    cat verkko-trio/6-rukki/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.paths.gaf | sed s/MAT/HAPLOTYPE1/g | sed s/PAT/HAPLOTYPE2/g > $flav/8-hicPipeline/rukki.paths.gaf && \\\n";
    print CMD "    cat verkko-trio/6-rukki/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.paths.tsv | sed s/MAT/HAPLOTYPE1/g | sed s/PAT/HAPLOTYPE2/g > $flav/8-hicPipeline/rukki.paths.tsv\n";
    print CMD "  fi\n";
    print CMD "\n";
    print CMD "  if [ ! -e $flav/assembly.fasta -a -e $flav/8-hicPipeline/rukki.paths.tsv ] ; then\n";
    print CMD "    $root/software/verkko/bin/verkko --graphaligner $root/software/graphaligner/bin/GraphAligner --slurm -d $flav \\\n";
    print CMD "      --ovb-run 8 32 32 \\\n";
    print CMD "      --screen human \\\n";
    #rint CMD "      --rdna-scaff \\\n";   #  Verkko 2.0 needs this.
    print CMD "      --hifi $hifi \\\n";
    print CMD "      --nano $nano \\\n";
    print CMD "      --hic1 $hic1 \\\n";
    print CMD "      --hic2 $hic2 \\\n";
    print CMD "    && \\\n";
    print CMD "    tar -C $flav               -cf - .snakemake    | gzip -9 > $flav/snakemake.tar.gz                   && rm -rf $flav/.snakemake && \\\n";
    print CMD "    tar -C $flav/8-hicPipeline -cf - .snakemake    | gzip -9 > $flav/8-hicPipeline/snakemake.tar.gz     && rm -rf $flav/8-hicPipeline/.snakemake && \\\n";
    print CMD "    tar -C $flav               -cf - batch-scripts | gzip -9 > $flav/batch-scripts.tar.gz               && rm -rf $flav/batch-scripts && \\\n";
    print CMD "    tar -C $flav/8-hicPipeline -cf - batch-scripts | gzip -9 > $flav/8-hicPipeline/batch-scripts.tar.gz && rm -rf $flav/8-hicPipeline/batch-scripts\n";
    print CMD "  fi\n";
    print CMD "fi\n";
    print CMD "\n";
    print CMD "rm -f $sdir/$flav.jid\n";
  }

  return "$sdir/$flav";
}

1;
