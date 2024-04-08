
###############################################################################
#
#  This file is part of a pipeline to (help) run Verkko assemblies on
#  numerous HPRC samples.
#
#  This is a 'United States Government Work', and is released in the public
#  domain.
#
##

package hprc::readFilter;
require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(filterReads);

use strict;
use warnings "all";
no  warnings "uninitialized";

use File::Basename;

use hprc::samples;
use hprc::assemble;


sub filterHiFi ($$$$) {
  my $file   = shift @_;   #  Full path to input
  my $baseIN = shift @_;   #  Full path without extensions
  my $name   = shift @_;   #  Just the name of the reads, no prefix, no path, no extensions
  my $submit = shift @_;   #    (the above is only for logging)

  my $bdir = dirname($baseIN) . "/hifi-cutadapt";   #  Throw output into hifi-cutadapt directory,
  my $bnam = basename($baseIN);                     #  using the original filename.

  system("mkdir -p $bdir");

  if (-e "$bdir/$bnam.jid") {
    print STDERR " - $name - RUNNING.\n";
    return;
  }
  if (-e "$bdir/$bnam.err") {
    print STDERR " - $name - CRASHED.\n";
    return;
  }
  if (-e "$bdir/$bnam.fasta.gz") {
    print STDERR " - $name - FINISHED.\n";
    return;
  }

  if (! -e "$bdir/$bnam.sh") {
    open(CMD, "> $bdir/$bnam.sh") or die "Failed to open '$bdir/$bnam.sh' for writing: $!";
    print CMD "#!/bin/sh\n";
    print CMD "#\n";
    print CMD "#SBATCH --cpus-per-task=24\n";
    print CMD "#SBATCH --mem=32g\n";
    print CMD "#SBATCH --time=4:00:00\n";
    print CMD "#SBATCH --output=$bdir/$bnam.err\n";
    print CMD "#SBATCH --job-name=cutada$name\n";
    print CMD "#\n";
    print CMD "\n";
    print CMD "set -o pipefail\n";
    print CMD "set -x\n";
    print CMD "\n";
    print CMD "cd $bdir\n";
    print CMD "\n";
    print CMD "export REF_CACHE=$root/samtools-ref-cache\n";
    print CMD "\n";
    print CMD "module load samtools\n";
    print CMD "module load cutadapt\n";
    print CMD "\n";
    print CMD "cutCPUs=2\n";
    print CMD "zipCPUs=\$(expr \$SLURM_JOB_CPUS_PER_NODE - 2)\n";
    print CMD "zipOpt='-9 -b 131072'\n";
    print CMD "\n";
    print CMD "if [ x\$zipCPUs = x ] ; then\n";
    print CMD "  echo 'WARNING: SLURM_JOB_CPUS_PER_NODE not set, using minimal CPUs instead.'\n";
    print CMD "  cutCPUs=1\n";
    print CMD "  zipCPUs=1\n";
    print CMD "  zipOpt='-1'\n";
    print CMD "fi\n";
    print CMD "\n";
    print CMD "\n";
    print CMD "$root/software/seqrequester/build/bin/seqrequester extract -fasta \\\n";
    print CMD "  '$file' \\\n";
    print CMD "| \\\n";
    print CMD "cutadapt --fasta --discard --revcomp -j \$cutCPUs -e 0.05 \\\n";
    print CMD " -b 'AAAAAAAAAAAAAAAAAATTAACGGAGGAGGAGGA;min_overlap=35' \\\n";
    print CMD " -b 'ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT;min_overlap=45' \\\n";
    print CMD " --json '$bnam.log.json' \\\n";
    print CMD " - \\\n";
    print CMD "| \\\n";
    print CMD "pigz \$zipOpt -p \$zipCPUs \\\n";
    print CMD " > '$bnam.fasta.gz' \\\n";
    print CMD "\n";
    print CMD "if [ $? != 0 ] ; then\n";
    print CMD "  echo 'Failed!'\n";
    print CMD "  rm $bnam.fasta.gz\n";
    print CMD "else\n";
    print CMD "  rm $bnam.err\n";
    print CMD "fi\n";
    print CMD "\n";
    print CMD "rm $bnam.jid\n";
    close(CMD);
  }

  if ($submit) {
    print STDERR " - $name - SUBMITTING.\n";
    system("sbatch $bdir/$bnam.sh > $bdir/$bnam.jid")   if ($submit);
  }
  else {
    print STDERR " - $name - NOT STARTED.\n";
  }
}


sub filterReads ($$$) {
  my $samp   = shift @_;
  my $types  = shift @_;
  my $opts   = shift @_;
  my $submit = exists $$opts{"submit"};

  foreach my $type (sort keys %$types) {
    my @files = getFiles($samp, $type);

    foreach my $file (@files) {
      my $cmd;
      my $base = $file;
      my $name;

      $base =~ s/.(fasta.gz|fastq.gz|sam|bam|cram)$//;
      $name =  $base;
      $name =~ s/^.*--//;

      if ($type eq 'hifi')     { filterHiFi($file, $base, $name, $submit); }
      if ($type eq 'ont')      {}
      if ($type eq 'hic')      {}
      if ($type eq 'ilmn')     {}
      if ($type eq 'ilmn-mat') {}
      if ($type eq 'ilmn-pat') {}
    }  #  Over files in type
  }    #  Over types
}

1;
