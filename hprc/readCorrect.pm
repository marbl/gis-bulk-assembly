
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

sub getOutputs($) {
  my $odir       =  shift @_;
  my $onam       =  "";
  my $hasHybrid = (`$rsoft/hifiasm/hifiasm 2>&1` =~ /--hf/);

  if ($hasHybrid) {
     $odir .= "/ont-hifiasm-hybrid-correct";#  Path in local filesystem (or empty if doesn't exist).
	 $onam  = "r10-hifiasm-correct.ec.ont.fq.gz";
  } else {
     $odir .= "/ont-hifiasm-correct";
	 $onam  = "r10-hifiasm-correct.ec.fq.gz";
  }
  return ($odir, $hasHybrid, $onam);
}

sub getCorrectedFiles($) {
  my $samp   = shift @_;
  my $type   = 'ont-r10';
  my %files = getFileMap($samp, $type, 1);
  my $input = scalar(keys %files);
  my ($odir, $hasHybrid, $onam) = getOutputs("$data/$samp");
  my $exprun   = "$odir/r10-hifiasm-correct.input.fastq.gz";
  my $expected = "$odir/$onam";

  return ((-e $expected || -e $exprun) && $input >= 1) ? "$expected" : "";
}

sub correctONTR10 ($$$$) {
  my $samp       =  shift @_;                       #  Sample name, for logging.
  my $inputs_ref =  shift @_;                       #  Read inputs
  my $orig       =  shift @_;
  my $onam       =  "r10-hifiasm-correct";
  my $submit     =  shift @_;                       #
  my @inputs = @$inputs_ref;

  # figure out if we have older hifiasm or version supporting hybrid correction
  my ($odir, $hasHybrid, $ig) = getOutputs($orig);
  my $hifiMissing;
  my $hifi;
  my $extraArgs;
  my $coverage = 1;
  if ($hasHybrid) {
     $hifi = getDownloadedFiles($samp, "hifi");   #  Return raw form of hifi data.
     $hifiMissing = (($hifi eq "") && (numFiles($samp, "hifi") > 0) ||
                       scalar(split ' ', $hifi) < numFiles($samp, "hifi"));
    if ($hifiMissing) {
       print "ERROR: $samp: have hybrid hifiasm from $rsoft/hifiasm/hifiasm but HiFi data is not available, not doing hybrid correction!\n" if $hifiMissing;
    }
    my %readTypes = $hifiMissing ? ( "ont-r10" => 1 ) : ( "hifi" => 1, "ont-r10" => 1 );
    $coverage = $hifiMissing ? 0 : getCoverage($samp, \%readTypes);                            # force failure when we're missing HiFi, can edit this and just get coverage and it will work as non-hybrid version
  }

  if ($coverage > 0 && ! -e "$odir/$onam.sh") {
    system("mkdir -p $odir");

    open(CMD, "> $odir/$onam.sh") or die "Failed to open '$odir/$onam.sh' for writing: $!";
    print CMD "#!/bin/bash\n";
    print CMD "#\n";
    if ($hasHybrid) { # requires more memory/time
       print CMD "#SBATCH --cpus-per-task=96\n";
       # we randomly use the normal or the largemem queue
       # we can't submit to both but this way let us run more jobs
       if (int(rand(2)) == 0) {
          print CMD "#SBATCH --partition=norm\n";
       } else {
          print CMD "#SBATCH --partition=largemem\n";
       }
       print CMD "#SBATCH --mem=720g\n";
       print CMD "#SBATCH --time=72:00:00\n";
    } else {
       print CMD "#SBATCH --cpus-per-task=72\n";
       print CMD "#SBATCH --partition=norm\n";
       print CMD "#SBATCH --mem=350g\n";
       print CMD "#SBATCH --time=120:00:00\n";
    }
    print CMD "#SBATCH --output=$odir/$onam.err\n";
    print CMD "#SBATCH --job-name=hifia$samp\n";
    print CMD "#\n";
    print CMD "\n";
    print CMD "set -o pipefail\n";
    print CMD "set -e\n";
    print CMD "set -x\n";
    print CMD "\n";
    print CMD "cd $odir\n";
    print CMD "\n";
    print CMD "export REF_CACHE=$ENV{'REF_CACHE'}\n";
    print CMD "\n";
    print CMD "module load samtools\n";
    print CMD "module load seqtk\n";
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
    if ($hasHybrid && !$hifiMissing) {
       print CMD "   (\n";
       print CMD "   for f in $hifi; do\n";
       print CMD "     if [[ \$f == *.bam ]]; then\n";
       print CMD "        samtools fastq \"\$f\"\n";
       print CMD "     elif [[ \$f == *.fastq.gz || \$f == *.fq.gz ]]; then\n";
       print CMD "        zcat \"\$f\"\n";
       print CMD "     elif [[ \$f == *.fastq.bz2 || \$f == *.fq.bz2 ]]; then\n";
       print CMD "        bzcat \"\$f\"\n";
       print CMD "     elif [[ \$f == *.fastq || \$f == *.fq ]]; then\n";
       print CMD "        cat \"\$f\"\n";
       print CMD "     else\n";
       print CMD "        echo \"Warning: unrecognized format for hifiInput \$f\" >&2\n";
       print CMD "     fi\n";
       print CMD "   done\n";
       print CMD "   ) | bgzip -@ \$zipCPUs \$zipOpt -c -I $odir/hifi.input.fastq.gz.gzi - > $odir/hifi.input.fastq.gz \n";
       $extraArgs="--hf $odir/hifi.input.fastq.gz \\\n        --hom-cov $coverage"; # --chn-occ 5"
    }

    print CMD "   (\n";
    print CMD "   for f in " . join(' ', @inputs) . "; do\n";
    print CMD "     if [[ \$f == *.bam ]]; then\n";
    print CMD "        samtools fastq \"\$f\"\n";
    print CMD "     elif [[ \$f == *.fastq.gz || \$f == *.fq.gz ]]; then\n";
    print CMD "        zcat \"\$f\"\n";
    print CMD "     elif [[ \$f == *.fastq.bz2 || \$f == *.fq.bz2 ]]; then\n";
    print CMD "        bzcat \"\$f\"\n";
    print CMD "     elif [[ \$f == *.fastq || \$f == *.fq ]]; then\n";
    print CMD "        cat \"\$f\"\n";
    print CMD "     else\n";
    print CMD "        echo \"Warning: unrecognized format for hifiInput \$f\" >&2\n";
    print CMD "     fi\n";
    print CMD "   done\n";
    print CMD "   ) | seqtk seq -L 10000 - | bgzip -@ \$zipCPUs \$zipOpt -c -I $onam.input.fastq.gz.gzi - > $onam.input.fastq.gz \n";

    print CMD "$rsoft/hifiasm/hifiasm -e --write-ec --ont $onam.input.fastq.gz \\\n";
    print CMD "        -o $onam.WORKING -t \$cutCPUs $extraArgs && \\\n";
    print CMD "        bgzip -@ \$zipCPUs \$zipOpt $onam.WORKING.ec.fq && \\\n";
    print CMD "        samtools faidx $onam.WORKING.ec.fq.gz && \\\n";
    if ($hasHybrid) {
       if (!$hifiMissing) {
           print CMD "        cat $onam.WORKING.ec.fq.gz.fai | grep     \"^m\" | awk '{print \$1}' > $onam.hifi.ids && \\\n";
           print CMD "        seqtk subseq $onam.WORKING.ec.fq.gz $onam.hifi.ids | bgzip -@ \$zipCPUs \$zipOpt -c -I $onam.ec.hifi.fq.gz.gzi > $onam.ec.hifi.fq.gz && \\\n";
       } else {
           print CMD "        touch $onam.hifi.ids && \\\n";
       }
       print CMD "        cat $onam.WORKING.ec.fq.gz.fai | grep  -v \"^m\" | awk '{print \$1}' > $onam.ont.ids  && \\\n";
       print CMD "        seqtk subseq $onam.WORKING.ec.fq.gz $onam.ont.ids  | bgzip -@ \$zipCPUs \$zipOpt -c -I $onam.ec.ont.fq.gz.gzi  > $onam.ec.ont.fq.gz  && \\\n";
    }
    print CMD "        mv $onam.WORKING.ec.fq.gz $onam.ec.fq.gz && \\\n";
    print CMD "        mv $onam.WORKING.ec.fq.gz.gzi $onam.ec.fq.gz.gzi && \\\n";
    print CMD "        mv $onam.WORKING.ec.fq.gz.fai $onam.ec.fq.gz.fai && \\\n";
    print CMD "\n";
    print CMD "if [ \$? != 0 ] ; then\n";
    print CMD "  echo 'Failed!'\n";
    print CMD "  rm -f $onam.ec.fq.gz\n";
    print CMD "  rm -f $onam.ec.fq\n";
    print CMD "else\n";
    print CMD "  rm -f $odir/hifi.input.fastq.gz\n";
    print CMD "  rm -f $odir/hifi.input.fastq.gz.gzi\n";
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
  elsif ($coverage <= 0)                   { print "$samp - MISSING-SUMMARY-FOR-COVERAGE\n"; }
  elsif (! $submit)                        { print "$samp - READY-TO-COMPUTE\n"; }
  else                                     { print "$samp - SUBMITTED\n"; system("sbatch $odir/$onam.sh > $odir/$onam.jid"); }
}


sub correctReads ($$) {
  my $samp   = shift @_;
  my $type   = 'ont-r10';
  my $opts   = shift @_;
  my $submit = exists $$opts{"submit"};

  my %files  = getFileMap($samp, $type, 1);
  my $input  = scalar(keys %files);
  my $idir   = "$data/$samp";
  my ($odir, $hasHybrid, $onam) = getOutputs($idir);
  my @missing = grep { ! -e $_ } (values %files);
  if (-e "$odir/$onam")        { print "$samp - FINISHED\n"; }
  elsif ($input == 0)  {   print "$samp - NO APPROPRIATE READ TYPE FOUND - CHECK FOR '$type' DATA\n";}
  elsif (@missing)  {   foreach my $missing (@missing) { print "$samp - NOT-FETCHED      - $missing\n"; } }
  else              {   correctONTR10($samp, [values %files], $idir, $submit); }
}

1;
