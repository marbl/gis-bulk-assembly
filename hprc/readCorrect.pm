
###############################################################################
#
#  This file is part of a pipeline to (help) run Verkko assemblies on
#  numerous HPRC samples.
#
#  This is a 'United States Government Work', and is released in the public
#  domain.
#
##

package hprc::readCorrect;
require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(correctReads getCorrectedFiles);

use strict;
use warnings "all";
no  warnings "uninitialized";

use File::Basename;

use hprc::samples;
use hprc::aws;

sub getCorrectedFiles($) {
  my $samp   = shift @_;
  my $type   = 'ont-r10';
  my %files = getFileMap($samp, $type, 1);
  my $input = scalar(keys %files);
  my $expected = "$data/$samp/ont-hifiasm-correct/r10-hifiasm-correct.ec.fq.gz";

  return (-e $expected && $input > 0) ? "$expected" : "";
}

sub correctONTR10 ($$$$) {
  my $samp       =  shift @_;                       #  Sample name, for logging.
  my $inputs_ref =  shift @_;                       #  Read inputs
  my $odir       =  (shift @_) . "/ont-hifiasm-correct";#  Path in local filesystem (or empty if doesn't exist).
  my $onam       =  "r10-hifiasm-correct";
  my $submit     =  shift @_;                       #
  my @inputs = @$inputs_ref;

  if (! -e "$odir/$onam.sh") {
    system("mkdir -p $odir");

    open(CMD, "> $odir/$onam.sh") or die "Failed to open '$odir/$onam.sh' for writing: $!";
    print CMD "#!/bin/bash\n";
    print CMD "#\n";
    print CMD "#SBATCH --cpus-per-task=72\n";
    print CMD "#SBATCH --partition=largemem\n";
    print CMD "#SBATCH --mem=500g\n";
    print CMD "#SBATCH --time=96:00:00\n";
    print CMD "#SBATCH --output=$odir/$onam.err\n";
    print CMD "#SBATCH --job-name=hifia$samp\n";
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
    print CMD "module load hifiasm\n";
    print CMD "\n";
    print CMD "cutCPUs=\$SLURM_JOB_CPUS_PER_NODE\n";
    print CMD "zipCPUs=\$(expr \$SLURM_JOB_CPUS_PER_NODE - 2)\n";
    print CMD "zipOpt='-l 9 -i'\n";
    print CMD "\n";
    print CMD "if [ x\$zipCPUs = x ] ; then\n";
    print CMD "  echo 'WARNING: SLURM_JOB_CPUS_PER_NODE not set, using minimal CPUs instead.'\n";
    print CMD "  cutCPUs=1\n";
    print CMD "  zipCPUs=1\n";
    print CMD "  zipOpt='-1'\n";
    print CMD "fi\n";
    print CMD "\n";
    print CMD "\n";
    print CMD "samtools merge -u - " . join(' ', @inputs) ." | samtools fastq - |seqtk seq -L 1000 - | bgzip -@ \$zipCPUs \$zipOpt -c -I $onam.input.fastq.gz.gzi - > $onam.input.fastq.gz && \\\n";
    print CMD "hifiasm -e --write-ec --ont -o $onam -t \$cutCPUs \\\n";
    print CMD "        $onam.input.fastq.gz && \\\n";
    print CMD "        bgzip -@ \$zipCPUs \$zipOpt $onam.ec.fq \n";
    print CMD "\n";
    print CMD "if [ \$? != 0 ] ; then\n";
    print CMD "  echo 'Failed!'\n";
    print CMD "  rm -f $onam.ec.fq.gz\n";
    print CMD "  rm -f $onam.ec.fq\n";
    print CMD "else\n";
    print CMD "  rm -f $onam.input.fastq.gz\n";
    print CMD "  rm -f $onam.input.fastq.gz.gzi\n";
    print CMD "  rm -f $onam.err\n";
    print CMD "fi\n";
    print CMD "\n";
    print CMD "rm -f $onam.*.bin\n";
    print CMD "rm -f $onam.jid\n";
    close(CMD);
  }

  my @missing = grep { ! -e $_ } @inputs;
  if    (-e "$odir/$onam.jid")             { print "$samp - RUNNING\n"; }
  elsif (-e "$odir/$onam.err")             { print "$samp - CRASHED\n"; }
  elsif (-e "$odir/$onam.ec.fq.gz")        { print "$samp - FINISHED\n"; }
  elsif (@missing)                         { foreach my $missing (@missing) { print "$samp - NOT-FETCHED      - $missing\n"; } }
  elsif (! $submit)                        { print "$samp - READY-TO-COMPUTE\n"; }
  else                                     { print "$samp - SUBMITTED\n"; system("sbatch $odir/$onam.sh > $odir/$onam.jid"); }
}


sub correctReads ($$) {
  my $samp   = shift @_;
  my $type   = 'ont-r10'; 
  my $opts   = shift @_;
  my $submit = exists $$opts{"submit"};

  my %files = getFileMap($samp, $type, 1);
  my $input = scalar(keys %files);
  my $idir = "$data/$samp/";
  my @missing = grep { ! -e $_ } (values %files);
  if ($input == 0)  {   print "$samp - NO APPROPRIATE READ TYPE FOUND - CHECK FOR '$type' DATA\n";}
  elsif (@missing)  {   foreach my $missing (@missing) { print "$samp - NOT-FETCHED      - $missing\n"; } }
  else              {   correctONTR10($samp, [values %files], $idir, $submit); }
}

1;
