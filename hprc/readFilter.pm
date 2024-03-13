
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

use hprc::samples;
use hprc::assemble;


sub filterHiFi ($$$$) {
  my $file   = shift @_;
  my $base   = shift @_;
  my $name   = shift @_;
  my $submit = shift @_;

  if (-e "$base.cutadapt.jid") {
    print STDERR " - $name - RUNNING.\n";
    return;
  }
  if (-e "$base.cutadapt.err") {
    print STDERR " - $name - CRASHED.\n";
    return;
  }
  if (-e "$base.cutadapt.fasta.gz") {
    print STDERR " - $name - FINISHED.\n";
    return;
  }

  if (! -e "$base.cutadapt.sh") {
    open(CMD, "> $base.cutadapt.sh") or die "Failed to open '$base.cutadapt.sh' for writing: $!";
    print CMD "#!/bin/sh\n";
    print CMD "#\n";
    print CMD "#SBATCH --cpus-per-task=2\n";
    print CMD "#SBATCH --mem=16g\n";
    print CMD "#SBATCH --time=4-0\n";
    print CMD "#SBATCH --output=$base.cutadapt.err\n";
    print CMD "#SBATCH --job-name=cutada$name\n";
    print CMD "#\n";
    print CMD "\n";
    print CMD "set -o pipefail\n";
    print CMD "set -x\n";
    print CMD "\n";
    print CMD "cd $root\n";
    print CMD "\n";
    print CMD "export REF_CACHE=$root/samtools-ref-cache\n";
    print CMD "\n";
    print CMD "module load samtools\n";
    print CMD "module load cutadapt\n";
    print CMD "\n";
    print CMD "cutCPUs=2\n";
    print CMD "zipCPUs=\$(expr \$SLURM_JOB_CPUS_PER_NODE - 2)\n";
    print CMD "\n";
    print CMD "$root/software/seqrequester/build/bin/seqrequester extract -fasta \\\n";
    print CMD "  '$file' \\\n";
    print CMD "| \\\n";
    print CMD "cutadapt --fasta --discard --revcomp -j \$cutCPUs -e 0.05 \\\n";
    print CMD " -b 'AAAAAAAAAAAAAAAAAATTAACGGAGGAGGAGGA;min_overlap=35' \\\n";
    print CMD " -b 'ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT;min_overlap=45' \\\n";
    print CMD " --json '$base.cutadapt.log.json' \\\n";
    print CMD " - \\\n";
    print CMD "| \\\n";
    print CMD "pigz -9 -b 131072 -p \$zipCPUs \\\n";
    print CMD " > '$base.cutadapt.fasta.gz' \\\n";
    print CMD "\n";
    print CMD "if [ $? != 0 ] ; then\n";
    print CMD "  echo 'Failed!'\n";
    print CMD "  rm $base.cutadapt.fasta.gz\n";
    print CMD "else\n";
    print CMD "  rm $base.cutadapt.err\n";
    print CMD "fi\n";
    print CMD "\n";
    print CMD "rm $base.cutadapt.jid\n";
    close(CMD);
  }

  print STDERR " - $name - SUBMITTING.\n";
  system("sbatch $base.cutadapt.sh > $base.cutadapt.jid")   if ($submit);
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

      $file =~ s/.bam/.fasta.xz/;

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
