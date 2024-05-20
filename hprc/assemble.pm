
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
@EXPORT = qw(isFinished cleanupAssembly archiveAssembly computeAssembly);

use strict;
use warnings "all";
no  warnings "uninitialized";

use hprc::samples;
use hprc::aws;

#
#  Returns 0 is the assembly is not currently running.
#   - no .jid file or slurm reports the job-id is not currently active.
#   - assembly could be complete or crashed.
#   - if .jid file exists, it will be removed.
#
#  Returns 1 if the assemble is definitely currently running.
#   - .jid file exists and slurm positively reports the job is active.
#
#  Returns 2 if the assembly is possibly running, by we can't tell for sure.
#   - every other case; manual intervention required.
#

sub isRunning ($$) {
  my $samp = shift @_;
  my $flav = shift @_;
  my $dirn = "$rasm/$samp`";

  my $jid;      #  Submitted job ID we're looking for.
  my $jobid;    #  Job ID we found.
  my $state;    #
  my $reason;   #
  my $out;      #

  my $ret = 0;                          #  No .jid file
  my $msg = "idle";                     #  0 - NOT RUNNING

  if (-e "$dirn/$flav.jid") {
    open(JID, "< $dirn/$flav.jid") or die "Failed to open '$dirn/$flav.jid' for reading: $!\n";
    my $jid = int(<JID>);
    close(JID);

    open(JOBS, "dashboard_cli jobs --noheader --tab --fields jobid,state,state_reason,std_out --jobid $jid|") or die "Failed to open 'dashboard_cli jobs' for reading: $!\n";
    while (<JOBS>) {
      ($jobid, $state, $reason, $out) = split '\s+', $_;

      last if ($jobid eq $jid);
    }
    close(JOBS);

    if     (!defined($jobid)) {         #  No report found for expected jobid
      ($ret, $msg) = (2, "orphan");     #  2 - POSSIBLY RUNNING
    }
    elsif (($state eq "COMPLETED") ||   #  Job completed successfully or not.
           ($state eq "CANCELLED") ||   #  0 - NOT RUNNING
           ($state eq "FAILED") ||
           ($state eq "TIMEOUT")) {
      ($ret, $msg) = (0, "failed");
    }
    elsif (($state eq "CONFIGURING") || #  Job is pending.
           ($state eq "PENDING")) {     #  1 - RUNNING
      ($ret, $msg) = (1, "pending");
    }
    elsif (($state eq "RUNNING") ||     #  Confirmed running (or at least under slurm control)
           ($state eq "COMPLETING") ||  #  1 - RUNNING
           ($state eq "STOPPED") ||
           ($state eq "SUSPENDED")) {
      ($ret, $msg) = (1, "running");
    }
    else {                              #  Unrecognized
      ($ret, $msg) = (2, "unknown");    #  2 - POSSIBLY RUNNING
    }
  }

  return wantarray() ? ($ret, $jobid, $state, $msg) : $ret;
}


#
#  Return 0 if the assembly outputs do not exist.
#  Return 1 if the assembly outputs do     exist; the assembly is complete.
#
sub isFinished ($$) {
  my $samp = shift @_;
  my $flav = shift @_;
  my $fini = 0;
  my $dirn = "$rasm/$samp";

  if     ($flav eq "canu-hifi") {
    $fini = 1   if (-e "$dirn/canu-hifi/asm.contigs.fasta");
  }
  elsif  ($flav eq "canu-trio") {
    $fini = 1   if ((-e "$dirn/canu-trio/asm-haplotypeMAT/asm.contigs.fasta") &&
                    (-e "$dirn/canu-trio/asm-haplotypePAT/asm.contigs.fasta"));
  }
  elsif (($flav eq "verkko-trio") ||
         ($flav eq "verkko-hi-c") ||
         ($flav eq "verkko-full")) {
    $fini = 1   if (-e "$dirn/$flav/assembly.fasta");
  }
  else {
    die "Unknown flavor $flav in isFinished().\n";
  }

  return $fini;
}

#
#  Return the size of a directory in GB.
#
sub getDirectorySize ($) {
  my $dd = shift @_;
  my $sm = 0;

  open(DU, "du -sm $dd |");
  while (<DU>) {
    $sm = $1   if (m/^\s*(\d+)\s+/);
  }
  close(DU);

  return $sm / 1024.0;
}

#
#  Submit functions.
#

sub createCanuTrio ($$$$$$) {
  my $samp  = shift @_;
  my $flav  = "canu-trio";
  my $nano  = shift @_;
  my $mati  = shift @_;
  my $pati  = shift @_;
  my $missi = shift @_;
  my $compl = shift @_;
  my $dirn  = "$rasm/$samp";

  if (!$missi && !$compl) {
    system("mkdir -p $dirn")  if (! -d $dirn);

    open(CMD, "> $dirn/$flav.sh") or die "Failed to open '$dirn/$flav.sh' for writing: $!\n";
    print CMD "#!/bin/sh\n";
    print CMD "\n";
    print CMD "cd $dirn\n";
    print CMD "\n";
    print CMD "export REF_CACHE=$ENV{'REF_CACHE'}\n";
    print CMD "\n";
    print CMD "module load samtools\n";
    print CMD "\n";
    print CMD "$root/software/canu/build/bin/canu -p asm -d verkko \\\n";
    print CMD "  genomeSize=3.1g useGrid=true gridOptionsJobName=ct$samp \\\n";
    print CMD "  gridOptionsOvl=\"-t 6-0\" \\\n";
    print CMD "  -haplotypeMAT $mati \\\n";
    print CMD "  -haplotypePAT $pati \\\n";
    print CMD "  -nanopore $nano \n";
    print CMD "\n";
    close (CMD);
  }

  return "$dirn/$flav";
}


sub createCanuHiFi ($$$$) {
  my $samp  = shift @_;
  my $flav  = "canu-hifi";
  my $hifi  = shift @_;
  my $missi = shift @_;
  my $compl = shift @_;
  my $dirn  = "$rasm/$samp";

  if (!$missi && !$compl) {
    system("mkdir -p $dirn")  if (! -d $dirn);

    open(CMD, "> $dirn/$flav.sh") or die "Failed to open '$dirn/$flav.sh' for writing: $!\n";
    print CMD "#!/bin/sh\n";
    print CMD "\n";
    print CMD "cd $dirn\n";
    print CMD "\n";
    print CMD "export REF_CACHE=$ENV{'REF_CACHE'}\n";
    print CMD "\n";
    print CMD "$root/software/canu/build/bin/canu -p asm -d verkko \\\n";
    print CMD "  genomeSize=3.1g useGrid=true gridOptionsJobName=ch$samp \\\n";
    print CMD "  gridOptionsOvl=\"-t 4-0\" \\\n";
    print CMD "  -pacbio-hifi $hifi \n";
    print CMD "\n";
    close (CMD);
  }

  return "$dirn/$flav";
}


sub createVerkkoBase ($$$$$) {
  my $samp  = shift @_;
  my $flav  = "verkko-trio";
  my $hifi  = shift @_;
  my $nano  = shift @_;
  my $missi = shift @_;
  my $compl = shift @_;
  my $dirn  = "$rasm/$samp";

  my ($mati, $pati) = locateHapmers($samp);

  if (!$missi && !$compl) {
    system("mkdir -p $dirn")  if (! -d $dirn);

    open(CMD, "> $dirn/$flav.sh") or die "Failed to open '$dirn/$flav.sh' for writing: $!\n";
    print CMD "#!/bin/sh\n";
    print CMD "#\n";
    print CMD "#SBATCH --cpus-per-task=2\n";
    print CMD "#SBATCH --mem=16g\n";
    print CMD "#SBATCH --time=4-0\n";
    print CMD "#SBATCH --output=$dirn/$flav.%j.log\n";
    print CMD "#SBATCH --job-name=va$samp\n";
    print CMD "#\n";
    print CMD "set -o pipefail\n";
    print CMD "set -x\n";
    print CMD "\n";
    print CMD "mkdir -p $dirn\n";
    print CMD "cd       $dirn\n";
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
    print CMD "#  Run initial trio-mode assembly in verkko/.  Move results to\n";
    print CMD "#  verkko-trio/ leaving verkko/ ready for a restart in Hi-C mode.\n";
    print CMD "#\n";
    print CMD "\n";
    print CMD "if [ ! -e $flav/assembly.fasta ] ; then\n";
    print CMD "  $root/software/verkko/bin/verkko --graphaligner $root/software/graphaligner/bin/GraphAligner --slurm -d verkko-base \\\n";
    print CMD "    --ovb-run 8 32 32 \\\n";
    print CMD "    --screen human \\\n";
    print CMD "    --hifi $hifi \\\n";
    print CMD "    --nano $nano \\\n";
    print CMD "    --hap-kmers $mati \\\n";
    print CMD "                $pati trio\\\n";
    print CMD "    5-untip\n";
    print CMD "\n";
    print CMD "  #  If successful, make one more output, cleanup and save the trio results.\n";
    print CMD "  if [ -e verkko/assembly.fasta ] ; then\n";
    print CMD "    #  We need contig fasta for the 'analysis' phase.  Make it now to avoid a race later.\n";
    print CMD "    if [ -e verkko/5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.gfa ] ; then\n";
    print CMD "      awk '/^S/ { print \">\"\$2; print \$3 }' \\\n";
    print CMD "        < verkko/5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.gfa \\\n";
    print CMD "      | fold -c \\\n";
    print CMD "        > verkko/5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.fasta\n";
    print CMD "    fi\n";
    print CMD "\n";
    print CMD "    rm -rf verkko/0-correction\n";
    print CMD "\n";
    print CMD "    mv verkko/6-layoutContigs $flav/\n";
    print CMD "    mv verkko/6-rukki         $flav/\n";
    print CMD "    mv verkko/7-consensus     $flav/\n";
    print CMD "    mv verkko/assembly*       $flav/\n";
    print CMD "    mv verkko/.snakemake      $flav/\n";
    print CMD "    mv verkko/batch-scripts   $flav/\n";
    print CMD "    mv verkko/snakemake.sh    $flav/\n";
    print CMD "    mv verkko/verkko.yml      $flav/\n";
    print CMD "\n";
    print CMD "    cd $flav/\n";
    print CMD "\n";
    print CMD "    tar -cf - .snakemake    | gzip -9 > snakemake.tar.gz     && rm -rf .snakemake\n";
    print CMD "    tar -cf - batch-scripts | gzip -9 > batch-scripts.tar.gz && rm -rf batch-scripts\n";
    print CMD "  fi\n";
    print CMD "fi\n";
    print CMD "\n";
    print CMD "rm -f $dirn/$flav.jid\n";
  }

  return "$dirn/$flav";
}


sub locateHapmers ($) {
  my $samp  = shift @_;
  my $mati  = "$root/hprc-data/$samp/hapmers/mati.hapmers.meryl";
  my $pati  = "$root/hprc-data/$samp/hapmers/pati.hapmers.meryl";
  my $missi = ((! -e "$mati/merylIndex") ||
               (! -e "$pati/merylIndex"));

  return wantarray() ? ($mati, $pati) : $missi;
}

sub createVerkkoTrio ($$$$$) {
  my $samp  = shift @_;
  my $flav  = "verkko-trio";
  my $hifi  = shift @_;
  my $nano  = shift @_;
  my $missi = shift @_;
  my $compl = shift @_;
  my $dirn  = "$rasm/$samp";

  my ($mati, $pati) = locateHapmers($samp);

  if (!$missi && !$compl) {
    system("mkdir -p $dirn")  if (! -d $dirn);

    open(CMD, "> $dirn/$flav.sh") or die "Failed to open '$dirn/$flav.sh' for writing: $!\n";
    print CMD "#!/bin/sh\n";
    print CMD "#\n";
    print CMD "#SBATCH --cpus-per-task=2\n";
    print CMD "#SBATCH --mem=16g\n";
    print CMD "#SBATCH --time=4-0\n";
    print CMD "#SBATCH --output=$dirn/$flav.%j.log\n";
    print CMD "#SBATCH --job-name=va$samp\n";
    print CMD "#\n";
    print CMD "set -o pipefail\n";
    print CMD "set -x\n";
    print CMD "\n";
    print CMD "mkdir -p $dirn\n";
    print CMD "cd       $dirn\n";
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
    print CMD "#  Run initial trio-mode assembly in verkko/.  Move results to\n";
    print CMD "#  verkko-trio/ leaving verkko/ ready for a restart in Hi-C mode.\n";
    print CMD "#\n";
    print CMD "\n";
    print CMD "if [ ! -e $flav/assembly.fasta ] ; then\n";
    print CMD "  $root/software/verkko/bin/verkko --graphaligner $root/software/graphaligner/bin/GraphAligner --slurm -d verkko-base \\\n";
    print CMD "    --ovb-run 8 32 32 \\\n";
    print CMD "    --screen human \\\n";
    print CMD "    --hifi $hifi \\\n";
    print CMD "    --nano $nano \\\n";
    print CMD "    --hap-kmers $mati \\\n";
    print CMD "                $pati trio\\\n";
    print CMD "    5-untip\n";
    print CMD "\n";
    print CMD "  #  If successful, make one more output, cleanup and save the trio results.\n";
    print CMD "  if [ -e verkko/assembly.fasta ] ; then\n";
    print CMD "    #  We need contig fasta for the 'analysis' phase.  Make it now to avoid a race later.\n";
    print CMD "    if [ -e verkko/5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.gfa ] ; then\n";
    print CMD "      awk '/^S/ { print \">\"\$2; print \$3 }' \\\n";
    print CMD "        < verkko/5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.gfa \\\n";
    print CMD "      | fold -c \\\n";
    print CMD "        > verkko/5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.fasta\n";
    print CMD "    fi\n";
    print CMD "\n";
    print CMD "    rm -rf verkko/0-correction\n";
    print CMD "\n";
    print CMD "    mv verkko/6-layoutContigs $flav/\n";
    print CMD "    mv verkko/6-rukki         $flav/\n";
    print CMD "    mv verkko/7-consensus     $flav/\n";
    print CMD "    mv verkko/assembly*       $flav/\n";
    print CMD "    mv verkko/.snakemake      $flav/\n";
    print CMD "    mv verkko/batch-scripts   $flav/\n";
    print CMD "    mv verkko/snakemake.sh    $flav/\n";
    print CMD "    mv verkko/verkko.yml      $flav/\n";
    print CMD "\n";
    print CMD "    cd $flav/\n";
    print CMD "\n";
    print CMD "    tar -cf - .snakemake    | gzip -9 > snakemake.tar.gz     && rm -rf .snakemake\n";
    print CMD "    tar -cf - batch-scripts | gzip -9 > batch-scripts.tar.gz && rm -rf batch-scripts\n";
    print CMD "  fi\n";
    print CMD "fi\n";
    print CMD "\n";
    print CMD "rm -f $dirn/$flav.jid\n";
  }

  return "$dirn/$flav";
}


sub createVerkkoHiC ($$$$$$$) {
  my $samp  = shift @_;
  my $flav  = "verkko-hi-c";
  my $hifi  = shift @_;
  my $nano  = shift @_;
  my $hic1  = shift @_;
  my $hic2  = shift @_;
  my $missi = shift @_;
  my $compl = shift @_;
  my $dirn  = "$rasm/$samp";

  if (!$missi && !$compl) {
    system("mkdir -p $dirn")  if (! -d $dirn);

    open(CMD, "> $dirn/$flav.sh") or die "Failed to open '$dirn/$flav.sh' for writing: $!\n";
    print CMD "#!/bin/sh\n";
    print CMD "#\n";
    print CMD "#SBATCH --cpus-per-task=2\n";
    print CMD "#SBATCH --mem=16g\n";
    print CMD "#SBATCH --time=4-0\n";
    print CMD "#SBATCH --output=$dirn/$flav.%j.log\n";
    print CMD "#SBATCH --job-name=vt$samp\n";
    print CMD "#\n";
    print CMD "set -o pipefail\n";
    print CMD "set -x\n";
    print CMD "\n";
    print CMD "mkdir -p $dirn\n";
    print CMD "cd       $dirn\n";
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
    print CMD "#  Given an initial trio-based assembly in verkko/ (with outputs moved to\n";
    print CMD "#  verkko-trio/) restart the assembly in Hi-C mode.  When that finishes,\n";
    print CMD "#  move files to verkko-hi-c/.\n";
    print CMD "#\n";
    print CMD "\n";
    print CMD "if [ ! -e $flav/assembly.fasta ] ; then\n";
    print CMD "  $root/software/verkko/bin/verkko --graphaligner $root/software/graphaligner/bin/GraphAligner --slurm -d $flav \\\n";
    print CMD "    --ovb-run 8 32 32 \\\n";
    print CMD "    --screen human \\\n";
    print CMD "    --hifi $hifi \\\n";
    print CMD "    --nano $nano \\\n";
    print CMD "    --hic1 $hic1 \\\n";
    print CMD "    --hic2 $hic2\n";
    print CMD "\n";
    print CMD "  mv verkko/6-layoutContigs $flav/\n";
    #rint CMD "  mv verkko/6-rukki         $flav/\n";
    print CMD "  mv verkko/7-consensus     $flav/\n";
    print CMD "  mv verkko/8-hicPipeline   $flav/\n";
    print CMD "  mv verkko/assembly*       $flav/\n";
    print CMD "  mv verkko/.snakemake      $flav/\n";
    print CMD "  mv verkko/batch-scripts   $flav/\n";
    print CMD "  mv verkko/snakemake.sh    $flav/\n";
    print CMD "  mv verkko/verkko.yml      $flav/\n";
    print CMD "\n";
    print CMD "  cd $flav/\n";
    print CMD "\n";
    print CMD "  tar -cf - .snakemake    | gzip -9 > snakemake.tar.gz     && rm -rf .snakemake\n";
    print CMD "  tar -cf - batch-scripts | gzip -9 > batch-scripts.tar.gz && rm -rf batch-scripts\n";
    print CMD "\n";
    print CMD "  cd 8-hicPipeline/\n";
    print CMD "\n";
    print CMD "  tar -cf - .snakemake    | gzip -9 > snakemake.tar.gz     && rm -rf .snakemake\n";
    print CMD "  tar -cf - batch-scripts | gzip -9 > batch-scripts.tar.gz && rm -rf batch-scripts\n";
    print CMD "fi\n";
    print CMD "\n";
    print CMD "rm -f $dirn/$flav.jid\n";
  }

  return "$dirn/$flav";
}



sub createVerkkoTrioHiC ($$$$$$$) {
  my $samp  = shift @_;
  my $flav  = "verkko-full";
  my $hifi  = shift @_;
  my $nano  = shift @_;
  my $hic1  = shift @_;
  my $hic2  = shift @_;
  my $missi = shift @_;
  my $compl = shift @_;
  my $dirn  = "$rasm/$samp";

  if (!$missi && !$compl) {
    system("mkdir -p $dirn")  if (! -d $dirn);

    open(CMD, "> $dirn/$flav.sh") or die "Failed to open '$dirn/$flav.sh' for writing: $!\n";
    print CMD "#!/bin/sh\n";
    print CMD "#\n";
    print CMD "#SBATCH --cpus-per-task=2\n";
    print CMD "#SBATCH --mem=16g\n";
    print CMD "#SBATCH --time=4-0\n";
    print CMD "#SBATCH --output=$dirn/$flav.%j.log\n";
    print CMD "#SBATCH --job-name=vh$samp\n";
    print CMD "#\n";
    print CMD "set -o pipefail\n";
    print CMD "set -x\n";
    print CMD "\n";
    print CMD "mkdir -p $dirn\n";
    print CMD "cd       $dirn\n";
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
    print CMD "#  Given an initial trio-based assembly in verkko/ (with outputs moved to\n";
    print CMD "#  verkko-trio/, and optionally a hi-c-based assembly in verkko-hi-c/) restart\n";
    print CMD "#  the assembly from the trio-mode path but using Hi-C for rDNA resolution.\n";
    print CMD "#  Outputs are left in verkko/.\n";
    print CMD "#\n";
    print CMD "\n";
    print CMD "if [ ! -e $flav/assembly.fasta ] ; then\n";
    print CMD "  if [ ! -e verkko/7-consensus ] ; then\n";
    print CMD "    mkdir verkko/7-consensus\n";
    print CMD "    ln -s ../../verkko-trio/7-consensus/ont_subset.fasta.gz  verkko/7-consensus/ont_subset.fasta.gz\n";
    print CMD "    ln -s ../../verkko-trio/7-consensus/ont_subset.id        verkko/7-consensus/ont_subset.id\n";
    print CMD "  fi\n";
    print CMD "\n";
    print CMD "  if [ ! -e verkko/8-hicPipeline ] ; then\n";
    print CMD "    mkdir verkko/8-hicPipeline\n";
    print CMD "    ln -s ../../verkko-trio/6-rukki/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.noseq.gfa  verkko/8-hicPipeline/unitigs.hpc.noseq.gfa\n";
    print CMD "    ln -s ../../verkko-trio/6-rukki/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.fasta      verkko/8-hicPipeline/unitigs.hpc.fasta\n";
    print CMD "  fi\n";
    print CMD "\n";
    print CMD "  if [ ! -e verkko/8-hicPipeline/rukki.paths.tsv ] ; then\n";
    print CMD "    $root/software/verkko/bin/verkko --graphaligner $root/software/graphaligner/bin/GraphAligner --slurm -d verkko --snakeopts '-U hicPhasing' \\\n";
    print CMD "      --ovb-run 8 32 32 \\\n";
    print CMD "      --screen human \\\n";
    #rint CMD "      --rdna-scaff \\\n";   #  Verkko 2.0 needs this.
    print CMD "      --hifi $hifi \\\n";
    print CMD "      --nano $nano \\\n";
    print CMD "      --hic1 $hic1 \\\n";
    print CMD "      --hic2 $hic2 \\\n";
    print CMD "    && \\\n";
    print CMD "    mv verkko/8-hicPipeline/hicverkko.colors.tsv                                                          verkko/8-hicPipeline/hicverkko.hiccolors.tsv \\\n";
    print CMD "    && \\\n";
    print CMD "    cp verkko-trio/6-rukki/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.colors.csv  verkko/8-hicPipeline/hicverkko.colors.tsv \\\n";
    print CMD "    && \\\n";
    print CMD "    rm -f verkko/8-hicPipeline/label* && \\\n";
    print CMD "    cp -p verkko-trio/6-rukki/label* verkko/8-hicPipeline/ \\\n";
    print CMD "    && \\\n";
    print CMD "    rm -f verkko/8-hicPipeline/rukki.paths.gaf && \\\n";
    print CMD "    rm -f verkko/8-hicPipeline/rukki.paths.tsv && \\\n";
    print CMD "    cat verkko-trio/6-rukki/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.paths.gaf | sed s/MAT/HAPLOTYPE1/g | sed s/PAT/HAPLOTYPE2/g > verkko/8-hicPipeline/rukki.paths.gaf && \\\n";
    print CMD "    cat verkko-trio/6-rukki/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.paths.tsv | sed s/MAT/HAPLOTYPE1/g | sed s/PAT/HAPLOTYPE2/g > verkko/8-hicPipeline/rukki.paths.tsv\n";
    print CMD "  fi\n";
    print CMD "\n";
    print CMD "  if [ ! -e verkko/assembly.fasta ] ; then\n";
    print CMD "    $root/software/verkko/bin/verkko --graphaligner $root/software/graphaligner/bin/GraphAligner --slurm -d verkko \\\n";
    print CMD "      --ovb-run 8 32 32 \\\n";
    print CMD "      --screen human \\\n";
    #rint CMD "      --rdna-scaff \\\n";   #  Verkko 2.0 needs this.
    print CMD "      --hifi $hifi \\\n";
    print CMD "      --nano $nano \\\n";
    print CMD "      --hic1 $hic1 \\\n";
    print CMD "      --hic2 $hic2\n";
    print CMD "  fi\n";
    print CMD "\n";
    print CMD "  if [ -e verkko/assembly.fasta ] ; then\n";
    print CMD "    mv verkko/6-layoutContigs $flav/\n";
    #rint CMD "    mv verkko/6-rukki         $flav/\n";
    print CMD "    mv verkko/7-consensus     $flav/\n";
    print CMD "    mv verkko/8-hicPipeline   $flav/\n";
    print CMD "    mv verkko/assembly*       $flav/\n";
    print CMD "    mv verkko/.snakemake      $flav/\n";
    print CMD "    mv verkko/batch-scripts   $flav/\n";
    print CMD "    mv verkko/snakemake.sh    $flav/\n";
    print CMD "    mv verkko/verkko.yml      $flav/\n";
    print CMD "\n";
    print CMD "    cd $flav\n";
    print CMD "    tar -cf - .snakemake    | gzip -9 > snakemake.tar.gz     && rm -rf .snakemake\n";
    print CMD "    tar -cf - batch-scripts | gzip -9 > batch-scripts.tar.gz && rm -rf batch-scripts\n";
    print CMD "\n";
    print CMD "    cd 8-hicPipeline/\n";
    print CMD "\n";
    print CMD "    tar -cf - .snakemake    | gzip -9 > snakemake.tar.gz     && rm -rf .snakemake\n";
    print CMD "    tar -cf - batch-scripts | gzip -9 > batch-scripts.tar.gz && rm -rf batch-scripts\n";
    print CMD "  fi\n";
    print CMD "fi\n";
    print CMD "\n";
    print CMD "rm -f $dirn/$flav.jid\n";
  }

  return "$dirn/$flav";
}


#
#  Remove intermediate files once assemblies that depend on them are complete.
#
#  In particular, the ONT 'split' directories and 0-correction cannot be removed
#  from verkko/ until ALL assemblies are complete.
#

sub cleanupAssembly ($$) {
  my $samp = shift @_;
  my $opts = shift @_;
  my @flavs;

  my $fintrio = isFinished($samp, "verkko-trio");
  my $finhic  = isFinished($samp, "verkko-hi-c");
  my $finfull = isFinished($samp, "verkko-full");
  my $finall  = $fintrio && $finhic && $finfull;

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
    my $dirn = "$rasm/$samp/$flav";

    if (! -e $dirn) {
      printf "%7s/%-11s - BEFORE: %6.1fGB - NOT STARTED.\n", $samp, $flav, 0.0;
      next;
    }

    printf "%7s/%-11s - BEFORE: ", $samp, $flav;   my $before = getDirectorySize($dirn);
    printf "%6.1fGB", $before;

    if ((($flav eq "verkko-trio") && (!$fintrio)) ||
        (($flav eq "verkko-hi-c") && (!$finhic)) ||
        (($flav eq "verkko-full") && (!$finfull)) ||
        (($flav eq "verkko")      && (!$finall))) {
      print " - INCOMPLETE.\n";
      next;
    }

    #  Save snakemake and script logs.  Each assembly command does this too.
    foreach my $d (qw(. 8-hicPipeline/final_contigs)) {
      if (-e "$dirn/$d/.snakemake")     { system("cd $dirn/$d && tar -cf snakemake-logs.tar .snakemake    && rm -rf .snakemake");    }
      if (-e "$dirn/$d/batch-scripts")  { system("cd $dirn/$d && tar -cf batch-scripts.tar  batch-scripts && rm -rf batch-scripts"); }
    }

    #  Remove correction intermediates and copies of ONT reads.
    if (-e "$dirn/0-correction")        { system("rm -rf $dirn/0-correction");        }
    if (-e "$dirn/3-align/split")       { system("rm -rf $dirn/3-align/split");       }
    if (-e "$dirn/3-alignTips/split")   { system("rm -rf $dirn/3-alignTips/split");   }
    if (-e "$dirn/8-hicPipeline/split") { system("rm -rf $dirn/8-hicPipeline/split"); }

    #  Remove consensus packages and job outputs.  Save logs.
    foreach my $d (qw(7-consensus 8-hicPipeline/final_contigs/7-consensus)) {
      if (-e "$dirn/$d/packages/part001.cnspack") { system("rm -rf $dirn/$d/packages/part*.cnspack"); }
      if (-e "$dirn/$d/packages/part001.fasta")   { system("rm -rf $dirn/$d/packages/part*.fasta");   }
      if (-e "$dirn/$d/packages")                 { system("cd $dirn/$d && tar -cf packages-logs.tar packages && rm -rf packages"); }
    }

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
    if (-e "$dirn/8-hicPipeline/align_bwa001.sh") { system("cd $dirn/8-hicPipeline && tar -cf align-bwa-logs.tar align_bwa* && rm -f align_bwa*"); }
    if (-e "$dirn/8-hicPipeline/unitigs.fasta.amb")  { system("rm -f $dirn/8-hicPipeline/unitigs.fasta.amb"); }
    if (-e "$dirn/8-hicPipeline/unitigs.fasta.ann")  { system("rm -f $dirn/8-hicPipeline/unitigs.fasta.ann"); }
    if (-e "$dirn/8-hicPipeline/unitigs.fasta.bwt")  { system("rm -f $dirn/8-hicPipeline/unitigs.fasta.bwt"); }
    if (-e "$dirn/8-hicPipeline/unitigs.fasta.pac")  { system("rm -f $dirn/8-hicPipeline/unitigs.fasta.pac"); }
    if (-e "$dirn/8-hicPipeline/unitigs.fasta.sa")  { system("rm -f $dirn/8-hicPipeline/unitigs.fasta.sa"); }

    printf " - AFTER: ";   my $after = getDirectorySize($dirn);
    printf "%6.1fGB (%.3f%%)\n", $after, 100.0 * $after / $before;
  }
}

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

        print STDOUT "  sbatch $path/$adir.tar.sh > $path/$adir.tar.jid\n";
        system("sbatch $path/$adir.tar.sh > $path/$adir.tar.jid 2>&1")   if ($subm);
      }
    }
  }
}

#
#  Main entry point
#

sub computeAssembly ($$) {
  my $samp    = shift @_;
  my $opts    = shift @_;
  my $flav    = $$opts{"flavor"};
  my $submit  = $$opts{"submit"};

  #  Check if output exists.

  my $compl = isFinished($samp, $flav);

  #  Check which inputs exist.

  my $hifi = getDownloadedFiles($samp, "hifi-cutadapt");   #  Return cutadapt form of hifi data.
  my $nano = getDownloadedFiles($samp, "ont");

  my $mati = getDownloadedFiles($samp, "mat-ilmn");
  my $pati = getDownloadedFiles($samp, "pat-ilmn");

  my $hic1 = getDownloadedFiles($samp, "hic1");
  my $hic2 = getDownloadedFiles($samp, "hic2");

  my $hifiMissing   =  (($hifi eq "") && (numFiles($samp, "hifi")     > 0));
  my $nanoMissing   =  (($nano eq "") && (numFiles($samp, "ont")      > 0));
  my $trioMissing   = ((($mati eq "") && (numFiles($samp, "mat-ilmn") > 0)) ||
                       (($pati eq "") && (numFiles($samp, "pat-ilmn") > 0)));
  my $hicMissing    = ((($hic1 eq "") && (numFiles($samp, "hic")      > 0)) ||
                       (($hic2 eq "") && (numFiles($samp, "hic")      > 0)));
  my $hapmerMissing = locateHapmers($samp);

  my $trioAsmMissing =   (! -e "$rasm/$samp/verkko-trio/assembly.fasta");

  #  Fail if all inputs are not created or downloaded.

  my $missing;

  if    ($flav eq "verkko-trio") {
    $missing .= "$samp/$flav -  can't start: missing hifi-cutadapt\n"     if ($hifiMissing);
    $missing .= "$samp/$flav -  can't start: missing ont\n"               if ($nanoMissing);
    $missing .= "$samp/$flav -  can't start: missing trio\n"              if ($trioMissing);
    $missing .= "$samp/$flav -  can't start: missing hapmer databases\n"  if ($hapmerMissing);
    #  NEED TO CHECK THAT TRIO DBs EXIST - done when launching
  }
  elsif ($flav eq "verkko-hi-c") {
    $missing .= "$samp/$flav -  can't start: missing hifi-cutadapt\n"     if ($hifiMissing);
    $missing .= "$samp/$flav -  can't start: missing ont\n"               if ($nanoMissing);
    $missing .= "$samp/$flav -  can't start: missing hic\n"               if ($hicMissing);
    $missing .= "$samp/$flav -  can't start: missing trio assembly\n"     if ($trioAsmMissing);
  }
  elsif ($flav eq "verkko-full") {
    $missing .= "$samp/$flav -  can't start: missing hifi-cutadapt\n"     if ($hifiMissing);
    $missing .= "$samp/$flav -  can't start: missing ont\n"               if ($nanoMissing);
    #missing .= "$samp/$flav -  can't start: missing trio\n"              if ($trioMissing);
    $missing .= "$samp/$flav -  can't start: missing hic\n"               if ($hicMissing);
    $missing .= "$samp/$flav -  can't start: missing trio assembly\n"     if ($trioAsmMissing);
  }
  elsif ($flav eq "canu-trio") {
    $missing .= "$samp/$flav -  can't start: missing ont\n"               if ($nanoMissing);
    $missing .= "$samp/$flav -  can't start: missing trio\n"              if ($trioMissing);
  }
  elsif ($flav eq "canu-hifi") {
    $missing .= "$samp/$flav -  can't start: missing hifi-cutadapt\n"     if ($hifiMissing);
  }
  elsif ($flav eq "") {
    $missing .= "$samp/{unknown} - can't start assembly without specifying a flavor.\n";
  }
  else {
    $missing .= "$samp/$flav - can't assemble unknown flavor '$flav'.\n";
  }

  #  Make a run script.

  my $scr;

  if ($flav eq "canu-trio")    { $scr = createCanuTrio     ($samp,        $nano, $mati, $pati, $missing, $compl); }
  if ($flav eq "canu-hifi")    { $scr = createCanuHiFi     ($samp, $hifi,                      $missing, $compl); }
  if ($flav eq "verkko-trio")  { $scr = createVerkkoTrio   ($samp, $hifi, $nano,               $missing, $compl); }   #  Trio DBs implicitly exist.
  if ($flav eq "verkko-hi-c")  { $scr = createVerkkoHiC    ($samp, $hifi, $nano, $hic1, $hic2, $missing, $compl); }
  if ($flav eq "verkko-full")  { $scr = createVerkkoTrioHiC($samp, $hifi, $nano, $hic1, $hic2, $missing, $compl); }

  if    ($compl)              { print "$samp/$flav - FINISHED\n";         }
  elsif (-e "$scr.jid")       { print "$samp/$flav - RUNNING\n";          }
  elsif (-e "$scr.err")       { print "$samp/$flav - CRASHED\n";          }
  elsif ($missing)            { print "$samp/$flav - MISSING-INPUTS\n";   }
  elsif (! $$opts{"submit"})  { print "$samp/$flav - READY-TO-COMPUTE\n"; }
  else                        { print "$samp/$flav - SUBMITTED\n"; system("sbatch $scr.sh > $scr.jid"); }
}

1;
