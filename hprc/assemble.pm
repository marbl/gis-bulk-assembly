
###############################################################################
#
#  This file is part of a pipeline to (help) run Verkko assemblies on
#  numerous HPRC samples.
#
#  This is a 'United States Government Work', and is released in the public
#  domain.
#
##

package hprc::assemble;
require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(numFiles getFiles startAssembly checkAssembly);

use strict;
use warnings "all";
no  warnings "uninitialized";

use hprc::samples;


sub numFiles ($$) {
  my $samp  = shift @_;
  my $type  = shift @_;
  my $fles  = $samples{$samp}{$type};
  return (defined($fles)) ? scalar(@$fles) : 0
}

sub getFiles ($$) {
  my $samp  = shift @_;
  my $type  = shift @_;
  my $hics;
  my $hici;
  my $subd;

  if      ($type eq "hic1") {     #  If asked for a specific Hi-C end, make
    $hics = "_R1_001.fastq.gz";   #  two extensions, one that we want to SAVE,
    $hici = "_R2_001.fastq.gz";   #  one that we want to INGORE.  We'll flag
    $type = "hic";                #  anything not in either of those as 'missing',
  } elsif ($type eq "hic2") {
    $hics = "_R2_001.fastq.gz";
    $hici = "_R1_001.fastq.gz";
    $type = "hic";
  }

  if ($type eq "hifi-cutadapt") {
    $subd = "hifi-cutadapt";
    $type = "hifi";
  }

  my $files = $samples{$samp}{$type};
  my @flist;
  my @missing;
  my @failed;

  foreach my $f (@$files) {
    my $awsf = $f;   #  So we don't accidentally corrupt the file list.
    my $locf = $f;
    my $locn = $f;

    $locf =~ s!/!--!g;
    $locf =~ s!^s3:----human-pangenomics--\w+--!$root/aws-data/$samp/!;

    $locn =~ s!^s3://human-pangenomics/\w+/!!;

    #  If a subdirectory is specified, insert it in the path AND use whatever
    #  extension exists in that directory.

    if (defined($subd)) {
      my $of = $locf;
      my $nf = $locf;

      $nf   =~ s/.(f(ast){0,1}[aq].gz|sam|bam|cram)$//i;  #  Replace existing extension
      $nf   =~ s!/$samp/!/$samp/$subd/!;                  #  Insert new subdirectory
      $locf = undef;

      foreach my $ext (qw(fasta.gz fa.gz fastq.gz fq.gz cram bam sam)) {
        if (-e "$nf.$ext") {
          $locf = "$nf.$ext";
          last;
        }
      }

      if (! -e $locf) {
        $locn =~ s/.(f(ast){0,1}[aq].gz|sam|bam|cram)$//i;   #  Make reported error name reflect
        $locn =~ s!/$samp/!/$samp/$subd/!;                   #  what we're searching for,
        
        push @missing, "$type\0$locn.(f(ast){0,1}[aq].gz|sam|bam|cram)";
      }
    }



    #  Some rather complicated rules to decide if the file is present, missing
    #  or if we're confused.

    if (!defined($hics) && (-e $locf)) {                            #  If not a Hi-C type,
      push @flist, $locf;                                           #  save if it exists.
    }                                                               #
    elsif (defined($hics) && (-e $locf) && ($locn =~ m/$hics$/)) {  #  If a Hi-C type, and we care
      push @flist, $locf;                                           #  about it, save if if exists.
    }                                                               #
    elsif (defined($hics)               && ($locn =~ m/$hici$/)) {  #  But ignore if it is the Hi-C
      ;                                                             #  end we don't care about.
    }                                                               #
    elsif (defined($hics)) {                                        #  Throw an error if it is a
      push @failed, "$type\0$locn";                                 #  file we don't recognize as
    }                                                               #  Hi-C.
    else {                                                          #
      push @missing, "$type\0$locn";                                #  Or declare we're missing the file.
    }
  }

  my $flist = join " \\\n         ", @flist;

  if (defined($hics) && (scalar(@failed) > 0)) {
    print STDERR "Can't launch the assembly, confused by Hi-C names:\n";

    foreach my $m (@failed) {
      my ($t, $f) = split '\0', $m;

      printf STDERR "  %8s - %s\n", $t, $f;
    }

    return undef;
  }
  else {
  }

  if (scalar(@missing) > 0) {
    print STDERR "Inputs not downloaded:\n";

    foreach my $m (@missing) {
      my ($t, $f) = split '\0', $m;

      printf STDERR "  %8s - %s\n", $t, $f;
    }

    return undef;
  }

  return wantarray() ? @flist : $flist;
}



sub isRunning ($$) {
  my $samp = shift @_;
  my $flav = shift @_;

  if (-e "$root/assemblies/$samp-$flav.jid") {
    my ($jobid, $state, $reason, $out) = (undef, undef, undef, undef);

    open(JID, "< $root/assemblies/$samp-$flav.jid") or die "Failed to open '$root/assemblies/$samp-$flav.jid' for reading: $!\n";
    my $jid = int(<JID>);
    close(JID);

    open(JOBS, "dashboard_cli jobs --noheader --tab --fields jobid,state,state_reason,std_out --jobid $jid|") or die "Failed to open 'dashboard_cli jobs' for reading: $!\n";
    while (<JOBS>) {
      ($jobid, $state, $reason, $out) = split '\s+', $_;

      last if ($jobid eq $jid);
    }
    close(JOBS);

    if (!defined($jobid)) {
      print STDERR "Didn't get a report for job $jid from Slurm.  Can't check if I can safely submit.\n";
      return 1;
    }

    if    (($state eq "CONFIGURING") ||
           ($state eq "PENDING")) {
      print STDERR "Job $jobid is in state $state.  Cannot submit it again.\n";
      return 1;
    }
    elsif (($state eq "RUNNING") ||
           ($state eq "COMPLETING") ||
           ($state eq "STOPPED") ||
           ($state eq "SUSPENDED")) {
      print STDERR "Job $jobid is in state $state.  Cannot submit it again.\n";
      return 1;
    }
    elsif (($state eq "COMPLETED") ||
           ($state eq "CANCELLED") ||
           ($state eq "FAILED") ||
           ($state eq "TIMEOUT")) {
      #print STDERR "Job $jobid is in state $state.  Will submit again.\n";
    }
    else {
      print STDERR "Job $jobid is in state $state.  Not sure if I can safely submit again.\n";
      return 1;
    }
  }

  return 0;
}

#
#  A rather boring test of python.  Might be useful when things break.
#

sub emitPythonTest () {
  print CMD "echo 'PYTHON VERSION:'\n";
  print CMD "command -v python\n";
  print CMD "python --version\n";
  print CMD "\n";
  print CMD "echo ''\n";
  print CMD "echo 'PYTHON LIBRARIES:'\n";
  print CMD "echo  > x.py 'from importlib import metadata'\n";
  print CMD "echo >> x.py 'for dist in metadata.distributions():'\n";
  print CMD "echo >> x.py '  print(dist._path)'\n";
  print CMD "python x.py\n";
  print CMD "\n";
  print CMD "echo ''\n";
  print CMD "echo 'NETWORKX VERSION'\n";
  print CMD "echo  > x.py 'from importlib.metadata import version'\n";
  print CMD "echo >> x.py 'version(\"networkx\")'\n";
  print CMD "python x.py\n";
  print CMD "\n";
  print CMD "echo ''\n";
  print CMD "echo 'NETWORKX TEST'\n";
  print CMD "echo  > x.py 'import networkx as nx'\n";
  print CMD "echo >> x.py 'G = nx.Graph()'\n";
  print CMD "echo >> x.py 'print(G)'\n";
  print CMD "python x.py\n";
  print CMD "\n";
  print CMD "rm x.py\n";
  print CMD "\n";
}

#
#  Helper functions - these just emit the verkko command; the
#  the boilerplate needed to actually run the command is in
#  the 'submit' functions later.
#

sub emitVerkkoTrio ($$$) {
  my $samp = shift @_;
  my $hifi = shift @_;
  my $nano = shift @_;

  my $mati = "$root/aws-data/$samp/hapmers/mati.hapmers.meryl";
  my $pati = "$root/aws-data/$samp/hapmers/pati.hapmers.meryl";

  if ((! -e "$mati/merylIndex") ||
      (! -e "$pati/merylIndex")) {
    print STDERR "Missing at least one hapmer database:\n";
    print STDERR "  $mati\n";
    print STDERR "  $pati\n";
    exit(1);
  }

  print CMD "\n";
  print CMD "\n";
  print CMD "#\n";
  print CMD "#  Run initial trio assembly, save results somewhere else,\n";
  print CMD "#  and setup for a Hi-C run.\n";
  print CMD "#\n";
  print CMD "\n";
  print CMD "if [ ! -e trio-outputs ] ; then\n";
  print CMD "  $root/software/verkko/bin/verkko --slurm -d . \\\n";
  print CMD "    --ovb-run 8 32 32 \\\n";
  print CMD "    --screen human \\\n";
  print CMD "    --hifi $hifi \\\n";
  print CMD "    --nano $nano \\\n";
  print CMD "    --hap-kmers $mati \\\n";
  print CMD "                $pati trio\n";
  print CMD "\n";
  print CMD "  mkdir trio-outputs\n";
  print CMD "\n";
  print CMD "  mv 6-layoutContigs trio-outputs/\n";
  print CMD "  mv 6-rukki         trio-outputs/\n";
  print CMD "  mv 7-consensus     trio-outputs/\n";
  print CMD "  mv assembly*       trio-outputs/\n";
  print CMD "  mv batch-scripts   trio-outputs/\n";
  print CMD "  mv snakemake.sh    trio-outputs/\n";
  print CMD "  mv verkko.yml      trio-outputs/\n";
  print CMD "fi\n";
}


sub emitVerkkoHiC ($$$$$) {
  my $samp = shift @_;
  my $hifi = shift @_;
  my $nano = shift @_;
  my $hic1 = shift @_;
  my $hic2 = shift @_;

  print CMD "\n";
  print CMD "\n";
  print CMD "#\n";
  print CMD "#  Run Hi-C\n";
  print CMD "#\n";
  print CMD "\n";
  print CMD "if [ ! -e trio-outputs/assembly.fasta ] ; then\n";
  print CMD "  echo \"Didn't find trio-assembly.\"\n";
  print CMD "  exit 1\n";
  print CMD "fi\n";
  print CMD "\n";
  print CMD "if [ ! -e 7-consensus ] ; then\n";
  print CMD "  mkdir 7-consensus\n";
  print CMD "  ln -s ../trio-outputs/7-consensus/ont_subset.fasta.gz  7-consensus/ont_subset.fasta.gz\n";
  print CMD "  ln -s ../trio-outputs/7-consensus/ont_subset.id        7-consensus/ont_subset.id\n";
  print CMD "fi\n";
  print CMD "\n";
  print CMD "if [ ! -e 8-hicPipeline ] ; then\n";
  print CMD "  mkdir 8-hicPipeline\n";
  print CMD "  ln -s ../trio-outputs/6-rukki/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.noseq.gfa  8-hicPipeline/unitigs.hpc.noseq.gfa\n";
  print CMD "  ln -s ../trio-outputs/6-rukki/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.fasta      8-hicPipeline/unitigs.hpc.fasta\n";
  print CMD "fi\n";
  print CMD "\n";
  print CMD "if [ ! -e assembly.fasta ] ; then\n";
  print CMD "  $root/software/verkko/bin/verkko --slurm -d . --snakeopts '-U hicPhasing' \\\n";
  print CMD "    --ovb-run 8 32 32 \\\n";
  print CMD "    --screen human \\\n";
  print CMD "    --rdna-scaff \\\n";
  print CMD "    --hifi $hifi \\\n";
  print CMD "    --nano $nano \\\n";
  print CMD "    --hic1 $hic1 \\\n";
  print CMD "    --hic2 $hic2\n";
  print CMD "\n";
  print CMD "  mv 8-hicPipeline/hicverkko.colors.tsv                                                                  8-hicPipeline/hicverkko.hiccolors.tsv\n";
  print CMD "  cp trio-outputs/6-rukki/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.colors.csv  8-hicPipeline/hicverkko.colors.tsv\n";
  print CMD "\n";
  print CMD "  cat trio-outputs/6-rukki/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.paths.gaf | sed s/MAT/HAPLOTYPE1/g | sed s/PAT/HAPLOTYPE2/g > 8-hicPipeline/rukki.paths.gaf\n";
  print CMD "  cat trio-outputs/6-rukki/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.paths.tsv | sed s/MAT/HAPLOTYPE1/g | sed s/PAT/HAPLOTYPE2/g > 8-hicPipeline/rukki.paths.tsv\n";
  print CMD "\n";
  print CMD "  cp -p trio-outputs/6-rukki/label* 8-hicPipeline/\n";
  print CMD "\n";
  print CMD "  $root/software/verkko/bin/verkko --slurm -d . \\\n";
  print CMD "    --ovb-run 8 32 32 \\\n";
  print CMD "    --screen human \\\n";
  print CMD "    --rdna-scaff \\\n";
  print CMD "    --hifi $hifi \\\n";
  print CMD "    --nano $nano \\\n";
  print CMD "    --hic1 $hic1 \\\n";
  print CMD "    --hic2 $hic2\n";
  print CMD "fi\n";
}

#
#  Submit functions.
#

sub submitCanuTrio ($$$$$) {
  my $samp = shift @_;
  my $nano = shift @_;
  my $mati = shift @_;
  my $pati = shift @_;
  my $subm = shift @_;

  open(CMD, "> $root/assemblies/$samp-canu-trio.sh") or die "Failed to open 'assemblies/$samp-canu-samp.sh' for writing: $!\n";
  print CMD "#!/bin/sh\n";
  print CMD "\n";
  print CMD "cd $root/assemblies/\n";
  print CMD "\n";
  print CMD "export REF_CACHE=$root/samtools-ref-cache\n";
  print CMD "\n";
  print CMD "module load samtools\n";
  print CMD "\n";
  print CMD "$root/software/canu/build/bin/canu -p asm -d $samp-canu-trio \\\n";
  print CMD "  genomeSize=3.1g useGrid=true gridOptionsJobName=ct$samp \\\n";
  print CMD "  gridOptionsOvl=\"-t 6-0\" \\\n";
  print CMD "  -haplotypeMAT $mati \\\n";
  print CMD "  -haplotypePAT $pati \\\n";
  print CMD "  -nanopore $nano \n";
  print CMD "\n";
  close (CMD);

  print STDOUT "sh $root/assemblies/$samp-canu-trio.sh > $root/assemblies/$samp-canu-trio.err\n";
  system("sh $root/assemblies/$samp-canu-trio.sh > $root/assemblies/$samp-canu-trio.err 2>&1")   if ($subm);
}


sub submitCanuHiFi ($$$) {
  my $samp = shift @_;
  my $hifi = shift @_;
  my $subm = shift @_;

  open(CMD, "> $root/assemblies/$samp-canu-hifi.sh") or die "Failed to open 'assemblies/$samp-canu-hifi.sh' for writing: $!\n";
  print CMD "#!/bin/sh\n";
  print CMD "\n";
  print CMD "cd $root/assemblies/\n";
  print CMD "\n";
  print CMD "export REF_CACHE=$root/samtools-ref-cache\n";
  print CMD "\n";
  print CMD "$root/software/canu/build/bin/canu -p asm -d $samp-canu-hifi \\\n";
  print CMD "  genomeSize=3.1g useGrid=true gridOptionsJobName=ch$samp \\\n";
  print CMD "  gridOptionsOvl=\"-t 4-0\" \\\n";
  print CMD "  -pacbio-hifi $hifi \n";
  print CMD "\n";
  close (CMD);

  print STDOUT "sh $root/assemblies/$samp-canu-hifi.sh > $root/assemblies/$samp-canu-hifi.err\n";
  system("sh $root/assemblies/$samp-canu-hifi.sh > $root/assemblies/$samp-canu-hifi.err 2>&1")   if ($subm);
}


sub submitVerkko ($$$$$$) {
  my $samp = shift @_;
  my $hifi = shift @_;
  my $nano = shift @_;
  my $hic1 = shift @_;
  my $hic2 = shift @_;
  my $subm = shift @_;

  open(CMD, "> $root/assemblies/$samp.sh") or die "Failed to open 'assemblies/$samp.sh' for writing: $!\n";
  print CMD "#!/bin/sh\n";
  print CMD "#\n";
  print CMD "#SBATCH --cpus-per-task=2\n";
  print CMD "#SBATCH --mem=16g\n";
  print CMD "#SBATCH --time=4-0\n";
  print CMD "#SBATCH --output=$root/assemblies/$samp.%j.log\n";
  print CMD "#SBATCH --job-name=va$samp\n";
  print CMD "#\n";
  print CMD "set -e\n";
  print CMD "set -x\n";
  print CMD "\n";
  print CMD "cd $root/assemblies/$samp\n";
  print CMD "\n";
  print CMD "export REF_CACHE=$root/samtools-ref-cache\n";
  print CMD "\n";
  print CMD "module load winnowmap\n";
  print CMD "module load mashmap\n";
  print CMD "module load samtools\n";
  print CMD "module load bwa\n";
  print CMD "module load seqtk\n";
  print CMD "\n";
  emitVerkkoTrio($samp, $hifi, $nano);
  emitVerkkoHiC ($samp, $hifi, $nano, $hic1, $hic2);

  print STDOUT "sbatch $root/assemblies/$samp.sh > $root/assemblies/$samp.jid\n";
  system("sbatch $root/assemblies/$samp.sh > $root/assemblies/$samp.jid 2>&1")   if ($subm);
}


sub submitVerkkoTrio ($$$$) {
  my $samp = shift @_;
  my $hifi = shift @_;
  my $nano = shift @_;
  my $subm = shift @_;

  open(CMD, "> $root/assemblies/$samp-trio.sh") or die "Failed to open 'assemblies/$samp-trio.sh' for writing: $!\n";
  print CMD "#!/bin/sh\n";
  print CMD "#\n";
  print CMD "#SBATCH --cpus-per-task=2\n";
  print CMD "#SBATCH --mem=16g\n";
  print CMD "#SBATCH --time=4-0\n";
  print CMD "#SBATCH --output=$root/assemblies/$samp-trio.%j.log\n";
  print CMD "#SBATCH --job-name=vt$samp\n";
  print CMD "#\n";
  print CMD "set -e\n";
  print CMD "set -x\n";
  print CMD "\n";
  print CMD "cd $root/assemblies/$samp\n";
  print CMD "\n";
  print CMD "export REF_CACHE=$root/samtools-ref-cache\n";
  print CMD "\n";
  print CMD "module load winnowmap\n";
  print CMD "module load mashmap\n";
  print CMD "module load samtools\n";
  print CMD "module load bwa\n";
  print CMD "module load seqtk\n";
  print CMD "\n";
  emitVerkkoTrio($samp, $hifi, $nano);
  #mitVerkkoHiC($hifi, $nano, $hic1, $hic2);

  print STDOUT "sbatch $root/assemblies/$samp-trio.sh > $root/assemblies/$samp-trio.jid\n";
  system("sbatch $root/assemblies/$samp-trio.sh > $root/assemblies/$samp-trio.jid 2>&1")   if ($subm);
}


sub submitVerkkoHiC ($$$$$$) {
  my $samp = shift @_;
  my $hifi = shift @_;
  my $nano = shift @_;
  my $hic1 = shift @_;
  my $hic2 = shift @_;
  my $subm = shift @_;

  open(CMD, "> $root/assemblies/$samp-hic.sh") or die "Failed to open 'assemblies/$samp-hic.sh' for writing: $!\n";
  print CMD "#!/bin/sh\n";
  print CMD "#\n";
  print CMD "#SBATCH --cpus-per-task=2\n";
  print CMD "#SBATCH --mem=16g\n";
  print CMD "#SBATCH --time=4-0\n";
  print CMD "#SBATCH --output=$root/assemblies/$samp-hic.%j.log\n";
  print CMD "#SBATCH --job-name=vh$samp\n";
  print CMD "#\n";
  print CMD "set -e\n";
  print CMD "set -x\n";
  print CMD "\n";
  print CMD "cd $root/assemblies/$samp\n";
  print CMD "\n";
  print CMD "export REF_CACHE=$root/samtools-ref-cache\n";
  print CMD "\n";
  print CMD "module load winnowmap\n";
  print CMD "module load mashmap\n";
  print CMD "module load samtools\n";
  print CMD "module load bwa\n";
  print CMD "module load seqtk\n";
  print CMD "\n";
  #mitVerkkoTrio($same, $hifi, $nano);
  emitVerkkoHiC ($samp, $hifi, $nano, $hic1, $hic2);

  print STDOUT "sbatch $root/assemblies/$samp-hic.sh > $root/assemblies/$samp-hic.jid\n";
  system("sbatch $root/assemblies/$samp-hic.sh > $root/assemblies/$samp-hic.jid 2>&1")   if ($subm);
}

#
#  Main entry point
#

sub startAssembly ($$) {
  my $samp   = shift @_;
  my $opts   = shift @_;
  my $submit = exists $$opts{"submit"};

  #  Make a place to work.

  if (! -d "$root/assemblies/$samp") {
    system("mkdir -p $root/assemblies/$samp");
  }

  #  Check that inputs exist.

  my $hifi = getFiles($samp, "hifi-cutadapt");
  my $nano = getFiles($samp, "ont");

  my $mati = getFiles($samp, "mat-ilmn");
  my $pati = getFiles($samp, "pat-ilmn");

  my $hic1 = getFiles($samp, "hic1");
  my $hic2 = getFiles($samp, "hic2");

  if (($hifi eq "") ||
      ($nano eq "") ||
      ($mati eq "") ||
      ($pati eq "") ||
      ($hic1 eq "") ||
      ($hic2 eq "")) {
    print STDERR "Inputs missing for sample $samp (check 'hprc.pl list --sample $samp').\n";
    return;
  }

  #  Decide which flavor of assembly to run.

  my @flavors;

  if (exists($$opts{"verkko-trio"}) &&
      exists($$opts{"verkko-hi-c"}))  {
    delete $$opts{"verkko-trio"};     #  Hi-C needs to run after trio, and they're
    delete $$opts{"verkko-hi-c"};     #  thus put into the same single run.
    push @flavors, "verkko";
  }

  if (exists($$opts{"canu-trio"}))    {  push @flavors, "canu-trio";    }
  if (exists($$opts{"canu-hifi"}))    {  push @flavors, "canu-hifi";    }
  if (exists($$opts{"verkko-trio"}))  {  push @flavors, "verkko-trio";  }
  if (exists($$opts{"verkko-hi-c"}))  {  push @flavors, "verkko-hi-c";  }

  if (scalar(@flavors) == 0)          {  push @flavors, "verkko";       }

  #  Make sure it isn't running.
  #    Does NOT work for Canu assemblies, or when a default verkko trio+hi-c
  #    assembly is started, then a verkko-trio is run.  So really of use only
  #    for the expected usual case of no assembly flavor specified on the
  #    command line.
  #
  foreach my $flav (@flavors) {
    if (isRunning($samp, $flav)) {
      return;
    }
  }

  #  Clean up a finished assembly?

  if ($$opts{'cleanup'}) {
    my $before = 0;
    my $after  = 0;

    if (! -e "assemblies/$samp/assembly.fasta") {
      print STDERR "ERROR:  Cannot clean up $samp; no assembly.fasta.\n";
      return;
    }

    printf STDERR "assemblies/$samp - BEFORE: ";

    open(DU, "du -sm assemblies/$samp |");
    while (<DU>) {
      $before = $1   if (m/^\s*(\d+)\s+/);
    }
    close(DU);

    printf STDERR "%.2fGB", $before / 1024.0;

    foreach my $d (qw(. trio-outputs 8-hicPipeline/final_contigs)) {
      if (-e "assemblies/$samp/$d/.snakemake") {
        system("cd assemblies/$samp/$d && tar -cf snakemake-logs.tar .snakemake && rm -rf .snakemake");
      }
      if (-e "assemblies/$samp/$d/batch-scripts") {
        system("cd assemblies/$samp/$d && tar -cf batch-scripts.tar batch-scripts && rm -rf batch-scripts");
      }
    }

    system("rm -rf assemblies/$samp/0-correction");

    foreach my $d (qw(3-align 3-alignTips 8-hicPipeline)) {
      system("rm -rf assemblies/$samp/$d/split");
    }

    foreach my $d (qw(7-consensus trio-outputs/7-consensus 8-hicPipeline/final_contigs/7-consensus)) {
      system("rm -rf assemblies/$samp/$d/packages/*.cnspack");
      system("rm -rf assemblies/$samp/$d/packages/*fasta");

      if (-e "assemblies/$samp/$d/packages") {
        system("cd assemblies/$samp/$d && tar -cf packages-logs.tar packages && rm -rf packages");
      }
    }

    if (-e "assemblies/$samp/8-hicPipeline/align_bwa001.sh") {
      system("cd assemblies/$samp/8-hicPipeline && tar -cf align-bwa-logs.tar align_bwa* && rm -f align_bwa*");
    }

    system("rm -f assemblies/$samp/8-hicPipeline/mapped???.bam");

    printf STDERR " - AFTER: ";

    open(DU, "du -sm assemblies/$samp |");
    while (<DU>) {
      $after = $1   if (m/^\s*(\d+)\s+/);
    }
    close(DU);

    printf STDERR "%03fGB (%.3f%%)\n", $after / 1024.0, 100.0 * $after / $before;

    return;
 }
 
  #  Launch.

  foreach my $flav (@flavors) {
    if ($flav eq "canu-trio")    { submitCanuTrio  ($samp,        $nano, $mati, $pati, $submit); }
    if ($flav eq "canu-hifi")    { submitCanuHiFi  ($samp, $hifi,                      $submit); }
    if ($flav eq "verkko")       { submitVerkko    ($samp, $hifi, $nano, $hic1, $hic2, $submit); }
    if ($flav eq "verkko-trio")  { submitVerkkoTrio($samp, $hifi, $nano,               $submit); }   #  Trio DBs implicitly exist.
    if ($flav eq "verkko-hi-c")  { submitVerkkoHiC ($samp, $hifi, $nano, $hic1, $hic2, $submit); }
  }
}

1;
