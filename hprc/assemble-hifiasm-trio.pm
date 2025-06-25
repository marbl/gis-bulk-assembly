
###############################################################################
#
#  This file is part of a pipeline to (help) run Hifiasm assemblies on
#  numerous HPRC samples.
#
#  This is a 'United States Government Work', and is released in the public
#  domain.
#
##

sub locateYakmers ($) {     #  In array context, returns the expected locations of the two hapmer
  my $samp  = shift @_;     #  databases.  In scalar context, returns false/true if those exist.
  my $mati  = "$data/$samp/yakmers/mati.yak";
  my $pati  = "$data/$samp/yakmers/pati.yak";
  my $missi = ((! -e "$mati") ||
               (! -e "$pati"));

  return wantarray() ? ($mati, $pati) : $missi;
}

sub createHifiasmTrio ($$$$$$$) {
  my $samp  = shift @_;
  my $flav  = "hifiasm-trio";
  my $hifi  = shift @_;
  my $nano  = shift @_;
  my $missi = shift @_;
  my $unava = shift @_;
  my $compl = shift @_;
  my $params= shift @_;
  my $sdir  = "$rasm/$samp";

  # hifiasm wants UL and HiC data to be comma-separated (but not hifi)
  $nano =~ s/ /,/g;

  if (!$missi && !$compl && !$unava && !-e "$sdir/$flav.sh") {
    system("mkdir -p $sdir")  if (! -d $sdir);

    open(CMD, "> $sdir/$flav.sh") or die "Failed to open '$sdir/$flav.sh' for writing: $!\n";
    print CMD "#!/bin/sh\n";
    print CMD "#\n";
    print CMD "#SBATCH --cpus-per-task=48\n";
    print CMD "#SBATCH --partition=largemem\n";
    print CMD "#SBATCH --mem=500g\n";
    print CMD "#SBATCH --time=6-0\n";
    print CMD "#SBATCH --output=$sdir/$flav.%j.log\n";
    print CMD "#SBATCH --job-name=hifitrio$samp\n";
    print CMD "#\n";
    print CMD "set -o pipefail\n";
    print CMD "set -e\n";
    print CMD "set -x\n";
    print CMD "\n";
    print CMD "mkdir -p $sdir\n";
    print CMD "cd       $sdir\n";
    print CMD "mkdir -p $flav\n";
    print CMD "\n";
    print CMD "export REF_CACHE=$ENV{'REF_CACHE'}\n";
    print CMD "\n";
    print CMD "module load python   # ensure we load newer python as we need 3.8+ for hic\n";
    print CMD "module load hifiasm\n";
    print CMD "module load winnowmap\n";
    print CMD "module load mashmap\n";
    print CMD "module load samtools\n";
    print CMD "module load bwa\n";
    print CMD "module load seqtk\n";
    print CMD "\n";
    print CMD "#\n";
    print CMD "if [ ! -e \"$flav/hifiasm.complete\" ]; then\n";
    print CMD "  $rsoft/hifiasm/hifiasm -d $flav/assembly \\\n";
    print CMD "    --dual-scaf --telo-m CCCTAA \\\n";
    print CMD "    -t \$SLURM_CPUS_PER_TASK    \\\n";
    print CMD "         $hifi \\\n";
    print CMD "    --ul $nano \\\n";
    print CMD "    --ul-cut 50000 \\\n";
    print CMD "    -1 $data/$samp/yakmers/pati.yak \\\n";
    print CMD "    -2 $data/$samp/yakmers/mati.yak \\\n";
    print CMD "    $params      \\\n";
    print CMD "  && \\\n";
    print CMD "  cat $flav/assembly.dip.hap1.p_ctg.gfa     |awk '{if (match(\$1, \"^S\")) { print \">\"\$2; print \$3}}' | fold -c > $flav/assembly.haplotype1.fasta && \\\n";
    print CMD "  cat $flav/assembly.dip.hap2.p_ctg.gfa     |awk '{if (match(\$1, \"^S\")) { print \">\"\$2; print \$3}}' | fold -c > $flav/assembly.haplotype2.fasta && \\\n";
    print CMD "  cat $flav/assembly.dip.hap[12].p_ctg.gfa  |awk '{if (match(\$1, \"^S\")) { print \">\"\$2; print \$3}}' | fold -c > $flav/assembly.fasta && \\\n";
    print CMD "  cat $flav/assembly.dip.hap[12].p_ctg.gfa  |awk '{if (match(\$1, \"^S\")) { print \">utig4-\"\$2; print \$3}}' | seqtk hpc - | fold -c > $flav/contigs.fasta && \\\n";
    print CMD "  cat $flav/assembly.dip.hap[12].p_ctg.gfa  |awk '{if (match(\$1, \"^S\")) { print \"path \"\$2\" utig4-\"\$2 } }' > $flav/assembly.scfmap && \\\n";
    print CMD "  cat $flav/assembly.dip.hap[12].p_ctg.gfa  |awk '{if (match(\$1, \"^S\")) { print \$1\"\\tutig4-\"\$2\"\\t*\\tLN:i:\"length(\$3); } else if (match(\$1, \"^L\")) print \$1\"\\tutig4-\"\$2\"\\t\"\$3\"\\tutig4-\"\$4\"\\t\"\$5\"\\t\"\$6 }' > $flav/assembly.homopolymer-compressed.noseq.gfa && \\\n";
    print CMD "  cat $flav/assembly.dip.hap[12].p_ctg.gfa  |awk '{if (match(\$1, \"^S\")) { gsub(/A+/, \"A\", \$3); gsub(/C+/, \"C\", \$3); gsub(/G+/, \"G\", \$3); gsub(/T+/, \"T\", \$3);  print \$1\"\\tutig4-\"\$2\"\\t\"\$3\"\\tLN:i:\"length(\$3); } else if (match(\$1, \"^L\")) print \$1\"\\tutig4-\"\$2\"\\t\"\$3\"\\tutig4-\"\$4\"\\t\"\$5\"\\t\"\$6 }' > $flav/assembly.homopolymer-compressed.gfa && \\\n";
    print CMD "  cat $flav/assembly.hic.hap[12].p_ctg.gfa  |awk '{if (NR == 1) print \"node\\tmat\\tpat\\tmat:pat\\tcolor\"; if (match(\$1, \"^S\")) { if (match(\$2, \"h1tg\")) { C=\"#FF8888\"; M=100000; P=0; } else if (match(\$2, \"h2tg\")) { C=\"#8888FF\"; M=0; P=100000; } else { C=\"#88FF88\"; M=0; P=0; } print \"utig4-\"\$2\"\\t\"M\"\\t\"P\"\\t\"M\":\"P\"\\t\"C}}' > $flav/assembly.colors.csv && \\\n";
    print CMD "  touch $flav/hifiasm.complete\n";
    print CMD "fi\n";
    print CMD "\n";
    print CMD "rm -f $sdir/$flav.jid\n";
  }

  return "$sdir/$flav";
}

1;
