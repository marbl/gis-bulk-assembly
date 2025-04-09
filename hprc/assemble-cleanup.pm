
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
#  Remove intermediate files once assemblies that depend on them are complete.
#
#  In particular, the ONT 'split' directories and 0-correction cannot be removed
#  from verkko/ until ALL assemblies are complete.
#


sub getDirectorySize ($) {    #
  my $dd = shift @_;          #  Return the size of a directory in GB.
  my $sm = 0;                 #

  open(DU, "du -sm $dd |");
  while (<DU>) {
    $sm = $1   if (m/^\s*(\d+)\s+/);
  }
  close(DU);

  return $sm / 1024.0;
}


sub cleanupAssembly ($$) {
  my $samp = shift @_;
  my $opts = shift @_;
  my @flavs;

  my $ta = ((numFiles($samp, "mat-ilmn") > 0) && (numFiles($samp, "pat-ilmn") > 0));
  my $ha =  (numFiles($samp, "hic")      > 0);
  my $th = $ta && $ha;

  my $finbase = isFinished($samp, "verkko-base");          #  Base is finished if it is finished.
  my $fintrio = isFinished($samp, "verkko-trio") || !$ta;  #  Trio is finished if it is finished or not available.
  my $finhic  = isFinished($samp, "verkko-hi-c") || !$ha;  #  Hi-C ... same
  my $finthic = isFinished($samp, "verkko-thic") || !$th;

  if (($$opts{"flavor"} ne "") && ($$opts{"flavor"} ne "verkko-full")) {
    push @flavs, $$opts{"flavor"};
  }
  else {
    push @flavs, "verkko-base";
    push @flavs, "verkko-trio"  if ($ta);
    push @flavs, "verkko-hi-c"  if ($ha);
    push @flavs, "verkko-thic"  if ($ta && $ha);
  }

  foreach my $flav (@flavs) {
    my $dirn = "$rasm/$samp/$flav";

    if (! -e $dirn) {
      printf "%7s/%-11s - BEFORE: %6.1fGB - NOT STARTED.\n", $samp, $flav, 0.0;
      next;
    }

    printf "%7s/%-11s - BEFORE: ", $samp, $flav;   my $before = getDirectorySize($dirn);
    printf "%6.1fGB", $before;

    if ((($flav eq "verkko-base") && (!$finbase)) ||
        (($flav eq "verkko-trio") && (!$fintrio)) ||
        (($flav eq "verkko-hi-c") && (!$finhic)) ||
        (($flav eq "verkko-thic") && (!$finthic))) {
      print " - INCOMPLETE.\n";
      next;
    }

    my $unused = ($finbase && $fintrio && $finhic && $finthic);

    if ($flav eq "verkko-base")   { $unused = ($finbase && $fintrio && $finhic && $finthic); }
    if ($flav eq "verkko-trio")   { $unused = ($finbase && $fintrio &&            $finthic); }
    if ($flav eq "verkko-hi-c")   { $unused = ($finbase &&             $finhic            ); }
    if ($flav eq "verkko-thic")   { $unused = ($finbase && $fintrio &&            $finthic); }

    if (!$unused) {
      print " - COMPLETE but in-use.\n";
      next;
    }

    #print " - COMPLETE.\n";   #  For testing!
    #next;

    #  Save snakemake and script logs.  Each assembly command does this too.
    foreach my $d (qw(. 8-hicPipeline/final_contigs)) {
      if (-e "$dirn/$d/.snakemake")     { system("cd $dirn/$d && tar -cf snakemake-logs.tar .snakemake    && rm -rf .snakemake");    }
      if (-e "$dirn/$d/batch-scripts")  { system("cd $dirn/$d && tar -cf batch-scripts.tar  batch-scripts && rm -rf batch-scripts"); }
    }

    #  Remove correction intermediates and copies of ONT reads.
    if (-e "$dirn/0-correction")        { system("rm -rf $dirn/0-correction");           }
    if (-e "$dirn/3-align/split")       { system("rm -rf $dirn/3-align/split");          }
    if (-e "$dirn/3-align/split")       { system("rm -rf $dirn/3-align/split.finished"); }
    if (-e "$dirn/3-alignTips/split")   { system("rm -rf $dirn/3-alignTips/split");      }
    if (-e "$dirn/8-hicPipeline/split") { system("rm -rf $dirn/8-hicPipeline/split");    }

    #  Remove consensus packages and job outputs.  Save logs.
    foreach my $d (qw(7-consensus 8-hicPipeline/final_contigs/7-consensus)) {
      if (-e "$dirn/$d/packages/part001.cnspack") { system("rm -rf $dirn/$d/packages/part*.cnspack"); }
      if (-e "$dirn/$d/packages/part001.fasta")   { system("rm -rf $dirn/$d/packages/part*.fasta");   }
      if (-e "$dirn/$d/packages")                 { system("cd $dirn/$d && tar -cf packages-logs.tar packages && rm -rf packages"); }
    }

    #  zip large gfa outputs
    if (-e "$dirn/2-processGraph/combined-edges-uniques.gfa")       { system("pigz $dirn/2-processGraph/combined-edges-uniques.gfa"); }
    if (-e "$dirn/2-processGraph/fixed-hifi-resolved.gfa")  { system("pigz $dirn/2-processGraph/fixed-hifi-resolved.gfa"); }
    if (-e "$dirn/2-processGraph/gapped-hifi-resolved.gfa") { system("pigz $dirn/2-processGraph/gapped-hifi-resolved.gfa"); }
    if (-e "$dirn/2-processGraph/gapped-once-hifi-resolved.gfa")    { system("pigz $dirn/2-processGraph/gapped-once-hifi-resolved.gfa"); }
    if (-e "$dirn/2-processGraph/gapped-twice-hifi-resolved.gfa")   { system("pigz $dirn/2-processGraph/gapped-twice-hifi-resolved.gfa"); }
    if (-e "$dirn/2-processGraph/unrolled-hifi-resolved.gfa")       { system("pigz $dirn/2-processGraph/unrolled-hifi-resolved.gfa"); }
    if (-e "$dirn/4-processONT/chopped-unitig-unrolled-hifi-resolved.gfa")  { system("pigz $dirn/4-processONT/chopped-unitig-unrolled-hifi-resolved.gfa"); }
    if (-e "$dirn/4-processONT/combined-edges-uniques.gfa")                 { system("pigz $dirn/4-processONT/combined-edges-uniques.gfa"); }
    if (-e "$dirn/4-processONT/connected.gfa")                              { system("pigz $dirn/4-processONT/connected.gfa"); }
    if (-e "$dirn/4-processONT/gapped-graphaln.gfa")                        { system("pigz $dirn/4-processONT/gapped-graphaln.gfa"); }
    if (-e "$dirn/4-processONT/gapped-once-unitig-unrolled-hifi-resolved.gfa")      { system("pigz $dirn/4-processONT/gapped-once-unitig-unrolled-hifi-resolved.gfa"); }
    if (-e "$dirn/4-processONT/gapped-thrice-unitig-unrolled-hifi-resolved.gfa")    { system("pigz $dirn/4-processONT/gapped-thrice-unitig-unrolled-hifi-resolved.gfa"); }
    if (-e "$dirn/4-processONT/gapped-twice-unitig-unrolled-hifi-resolved.gfa")     { system("pigz $dirn/4-processONT/gapped-twice-unitig-unrolled-hifi-resolved.gfa"); }
    if (-e "$dirn/4-processONT/gapped-unitig-unrolled-hifi-resolved.gfa")   { system("pigz $dirn/4-processONT/gapped-unitig-unrolled-hifi-resolved.gfa"); }
    if (-e "$dirn/4-processONT/normal-connected.gfa")       { system("pigz $dirn/4-processONT/normal-connected.gfa"); }
    if (-e "$dirn/4-processONT/ont-resolved-graph.gfa")     { system("pigz $dirn/4-processONT/ont-resolved-graph.gfa"); }
    if (-e "$dirn/4-processONT/unrolled-ont-resolved.gfa")  { system("pigz $dirn/4-processONT/unrolled-ont-resolved.gfa"); }
    if (-e "$dirn/5-untip/combined-edges-1.gfa")            { system("pigz $dirn/5-untip/combined-edges-1.gfa"); }
    if (-e "$dirn/5-untip/combined-edges-final.gfa")        { system("pigz $dirn/5-untip/combined-edges-final.gfa"); }
    if (-e "$dirn/5-untip/connected-tip.gfa")               { system("pigz $dirn/5-untip/connected-tip.gfa"); }
    if (-e "$dirn/5-untip/popped-connected-tip.gfa")        { system("pigz $dirn/5-untip/popped-connected-tip.gfa"); }
    if (-e "$dirn/5-untip/popped-unitig-unrolled-popped-connected-tip.gfa") { system("pigz $dirn/5-untip/popped-unitig-unrolled-popped-connected-tip.gfa"); }
    if (-e "$dirn/5-untip/unitig-connected-tip.gfa")        { system("pigz $dirn/5-untip/unitig-connected-tip.gfa"); }
    if (-e "$dirn/5-untip/unitig-unrolled-popped-connected-tip.gfa")        { system("pigz $dirn/5-untip/unitig-unrolled-popped-connected-tip.gfa"); }
    if (-e "$dirn/5-untip/unrolled-popped-connected-tip.gfa")       { system("pigz $dirn/5-untip/unrolled-popped-connected-tip.gfa"); }
    if (-e "$dirn/5-untip/unrolled-popped-unitig-unrolled-popped-connected-tip.gfa")        { system("pigz $dirn/5-untip/unrolled-popped-unitig-unrolled-popped-connected-tip.gfa"); }
    if (-e "$dirn/4-processONT/alns-cut.gaf")                                               { system("pigz $dirn/4-processONT/alns-cut.gaf");}
    if (-e "$dirn/4-processONT/alns-cut-trim.gaf")                                          { system("pigz $dirn/4-processONT/alns-cut-trim.gaf");}
    if (-e "$dirn/4-processONT/alns-graphalgn-nogap-1.gaf")                                 { system("pigz $dirn/4-processONT/alns-graphalgn-nogap-1.gaf");}
    if (-e "$dirn/4-processONT/alns-ont-filter-trim.gaf")                                   { system("pigz $dirn/4-processONT/alns-ont-filter-trim.gaf");}
    if (-e "$dirn/4-processONT/alns-ont-nogap-1.gaf")                                       { system("pigz $dirn/4-processONT/alns-ont-nogap-1.gaf");}
    if (-e "$dirn/4-processONT/alns-ont-nogap-2.gaf")                                       { system("pigz $dirn/4-processONT/alns-ont-nogap-2.gaf");}
    if (-e "$dirn/4-processONT/alns-ont-nogap-3.gaf")                                       { system("pigz $dirn/4-processONT/alns-ont-nogap-3.gaf");}
    if (-e "$dirn/4-processONT/alns-ont-nogap.gaf")                                         { system("pigz $dirn/4-processONT/alns-ont-nogap.gaf");}
    if (-e "$dirn/4-processONT/alns-trimmed.gaf")                                           { system("pigz $dirn/4-processONT/alns-trimmed.gaf");}
    if (-e "$dirn/4-processONT/fake-ont-alns.gaf")                                          { system("pigz $dirn/4-processONT/fake-ont-alns.gaf");}

    #  Remove consensus outputs that are copied to the assembly output directory.
    foreach my $f (qw(assembly.homopolymer-compressed.layout
                      assembly.disconnected.fasta
                      assembly.unassigned.fasta
                      assembly.fasta
                      assembly.haplotype1.fasta
                      assembly.haplotype2.fasta
                      assembly.ebv.exemplar.fasta  assembly.ebv.fasta
                      assembly.mito.exemplar.fasta assembly.mito.fasta
                      assembly.rdna.exemplar.fasta assembly.rdna.fasta
                      combined.fasta
                      unitig-popped.fasta
                      unitig-popped.haplotype1.fasta
                      unitig-popped.haplotype2.fasta
                      unitig-popped.unassigned.fasta)) {
      if (-e "$dirn/7-consensus/$f")                               { system("rm -f $dirn/7-consensus/$f"); }
      if (-e "$dirn/8-hicPipeline/final_contigs/$f")               { system("rm -f $dirn/8-hicPipeline/final_contigs/$f"); }
      if (-e "$dirn/8-hicPipeline/final_contigs/7-consensus/$f")   { system("rm -f $dirn/8-hicPipeline/final_contigs/7-consensus/$f"); }
    }

    #  Save bwa logs and scripts, but wipe the job outputs and databases.
    if (-e "$dirn/8-hicPipeline/mapped001.bam")   { system("rm -f $dirn/8-hicPipeline/mapped???.bam"); }
    if (-e "$dirn/8-hicPipeline/mapped001_nodefiltered.bam") { system("rm -f $dirn/8-hicPipeline/mapped???_nodefiltered.bam"); }
    if (-e "$dirn/8-hicPipeline/align_bwa001.sh") { system("cd $dirn/8-hicPipeline && tar -cf align-bwa-logs.tar align_bwa* && rm -f align_bwa*"); }
    if (-e "$dirn/8-hicPipeline/unitigs.fasta.amb")  { system("rm -f $dirn/8-hicPipeline/unitigs.fasta.amb"); }
    if (-e "$dirn/8-hicPipeline/unitigs.fasta.ann")  { system("rm -f $dirn/8-hicPipeline/unitigs.fasta.ann"); }
    if (-e "$dirn/8-hicPipeline/unitigs.fasta.bwt")  { system("rm -f $dirn/8-hicPipeline/unitigs.fasta.bwt"); }
    if (-e "$dirn/8-hicPipeline/unitigs.fasta.pac")  { system("rm -f $dirn/8-hicPipeline/unitigs.fasta.pac"); }
    if (-e "$dirn/8-hicPipeline/unitigs.fasta.sa")   { system("rm -f $dirn/8-hicPipeline/unitigs.fasta.sa"); }
    if (-e "$dirn/8-hicPipeline/unitigs.fasta")      { system("pigz $dirn/8-hicPipeline/unitigs.fasta"); }
    if (-e "$dirn/8-hicPipeline/paths.hpc.fasta")    { system("rm -f $dirn/8-hicPipeline/paths.hpc.fasta"); }

    printf " - AFTER: ";   my $after = getDirectorySize($dirn);
    printf "%6.1fGB (%.3f%%)\n", $after, 100.0 * ($before > 0 ? $after / $before : 1);
  }
}

1;
