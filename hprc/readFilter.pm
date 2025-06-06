
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
use hprc::aws;


sub filterHiFi ($$$$$$$) {
  my $samp   =  shift @_;                       #  Sample name, for logging.
  my $type   =  shift @_;                       #  Read type to process, for logging.
  my $locf   =  shift @_;                       #  Path in local filesystem (or empty if doesn't exist).
  my $odir   = (shift @_) . "/hifi-cutadapt";   #  Path to desired output directory.
  my $onam   =  shift @_;                       #  Desired output name, sans suffixes.
  my $lnam   =  shift @_;                       #  Name useful for logging.
  my $submit =  shift @_;                       #

  if ((-e $locf) && (! -e "$odir/$onam.sh")) {
    system("mkdir -p $odir");

    open(CMD, "> $odir/$onam.sh") or die "Failed to open '$odir/$onam.sh' for writing: $!";
    print CMD "#!/bin/sh\n";
    print CMD "#\n";
    print CMD "#SBATCH --cpus-per-task=24\n";
    print CMD "#SBATCH --mem=32g\n";
    print CMD "#SBATCH --time=4:00:00\n";
    print CMD "#SBATCH --output=$odir/$onam.err\n";
    print CMD "#SBATCH --job-name=cutada$lnam\n";
    print CMD "#\n";
    print CMD "\n";
    print CMD "set -o pipefail\n";
    print CMD "set -x\n";
    print CMD "\n";
    print CMD "cd $odir\n";
    print CMD "\n";
    print CMD "export REF_CACHE=$ENV{'REF_CACHE'}\n";
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
    print CMD "  '$locf' \\\n";
    print CMD "| \\\n";
    print CMD "cutadapt --fasta --discard --revcomp -j \$cutCPUs -e 0.05 \\\n";
    print CMD " -b 'AAAAAAAAAAAAAAAAAATTAACGGAGGAGGAGGA;min_overlap=35' \\\n";
    print CMD " -b 'ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT;min_overlap=45' \\\n";
    print CMD " --json '$onam.log.json' \\\n";
    print CMD " - \\\n";
    print CMD "| \\\n";
    print CMD "pigz \$zipOpt -p \$zipCPUs \\\n";
    print CMD " > '$onam.WORKING.fasta.gz' \\\n";
    print CMD "\n";
    print CMD "if [ $? != 0 ] ; then\n";
    print CMD "  echo 'Failed!'\n";
    print CMD "  rm $onam.WORKING.fasta.gz\n";
    print CMD "else\n";
    print CMD "  mv $onam.WORKING.fasta.gz $onam.fasta.gz\n";
    print CMD "  rm $onam.err\n";
    print CMD "fi\n";
    print CMD "\n";
    print CMD "rm $onam.jid\n";
    close(CMD);
  }

  if    (-e "$odir/$onam.jid")      { print "$samp/$type - RUNNING          - $lnam\n"; }
  elsif (-e "$odir/$onam.err")      { print "$samp/$type - CRASHED          - $lnam\n"; }
  elsif (-e "$odir/$onam.fasta.gz") { print "$samp/$type - FINISHED         - $lnam\n"; }
  elsif (! -e $locf)                { print "$samp/$type - NOT-FETCHED      - $lnam\n"; }
  elsif (! $submit)                 { print "$samp/$type - READY-TO-COMPUTE - $lnam\n"; }
  else                              { print "$samp/$type - SUBMITTED        - $lnam\n"; system("sbatch $odir/$onam.sh > $odir/$onam.jid"); }
}


sub filterReads ($$$) {
  my $samp   = shift @_;
  my $types  = shift @_;
  my $opts   = shift @_;
  my $submit = exists $$opts{"submit"};

  foreach my $type (sort keys %$types) {
    my %files = getFileMap($samp, $type, 1);

    foreach my $locf (values %files) {
      my $idir = dirname($locf);   #  Path to input file (converted to odir in filter function).
      my $bnam = basename($locf);  #  Name of input file.
      my $lnam = basename($locf);  #  Name for logging.

      $bnam =~ s/.(fasta.gz|fastq.gz|sam|bam|cram|fq.gz|fa.gz)$//;   #  Strip extensions from base name.

      $lnam =~ s!/!--!;                                  #  Change actual components to virtual,
      $lnam =~ s!^.*--!!;                                #  then erase all but the last component,
      $lnam =~ s!.(fasta.gz|fastq.gz|sam|bam|cram|fq.gz|fa.gz)$!!;   #  and finally erase the extension.

      if ($type eq 'hifi')     { filterHiFi($samp, $type, $locf, $idir, $bnam, $lnam, $submit); }
      if ($type eq 'ont')      {}
      if ($type eq 'hic')      {}
      if ($type eq 'ilmn')     {}
      if ($type eq 'ilmn-mat') {}
      if ($type eq 'ilmn-pat') {}
    }  #  Over files in type
  }    #  Over types
}

1;
