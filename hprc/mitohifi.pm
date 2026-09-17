
###############################################################################
#
#  This file is part of a pipeline to (help) run Verkko assemblies on
#  numerous HPRC samples.
#
#  This is a 'United States Government Work', and is released in the public
#  domain.
#
##

package hprc::mitohifi;
require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(startMitohifiAnalysis);

use strict;
use warnings "all";
no  warnings "uninitialized";

use hprc::aws;
use hprc::assemble;
use hprc::samples;

sub startMitohifiAnalysis($$) {
  my $samp   = shift @_;
  my $opts   = shift @_;
  my $flav   = $$opts{"flavor"};
  my $submit = $$opts{"submit"};
  my $diro   = "$rasm/$samp/$flav/mitohifi";

  #  Check if outputs exist.

  my $finished = -e "$diro/mitohifi.out";

  #  Check that inputs exist.

  my $ready     = isFinished($samp, $flav);
  my $hifi      = getDownloadedFiles($samp, "hifi-cutadapt");   #  Return cutadapt form of hifi data.
  my $unavail   = (($hifi eq "") && (numFiles($samp, "hifi-cutadapt") > 0) ||
                      scalar(split ' ', $hifi) < numFiles($samp, "hifi-cutadapt"));

  if ($ready && !$finished) {
    system("mkdir -p $diro")   if (! -d "$diro");

    open(CMD, "> $diro/mitohifi.sh") or die "Failed to open '$diro/mitohifi.sh' for writing: $!\n";
    print CMD "#!/bin/sh\n";
    print CMD "#\n";
    print CMD "#SBATCH --cpus-per-task=4\n";
    print CMD "#SBATCH --mem=50g\n";
    print CMD "#SBATCH --time=1:00:00\n";
    print CMD "#SBATCH --partition=quick\n";
    print CMD "#SBATCH --output=$diro/mitohifi.%j.err\n";
    print CMD "#SBATCH --job-name=mit$samp\n";
    print CMD "#\n";
    print CMD "set -e\n";
    print CMD "set -x\n";
    print CMD "\n";
    print CMD "mkdir -p $diro\n";
    print CMD "cd       $diro\n";
    print CMD "\n";
    print CMD "module load minimap2\n";
    print CMD "module load samtools\n";
    print CMD "module load seqtk\n";
    print CMD "\n";
    print CMD "\n";
    print CMD "if [ ! -e mitohifi.out ] ; then\n";
    # get the tig ID
    print CMD "TIG_ID=`head -n 1 ../assembly.mito.exemplar.fasta |awk '{print substr(\$1, 2, length(\$1))}'`\n";
    # extract the reads
    print CMD "python $rsoft/mitohifi/src/extract_reads_from_verkko.py \$TIG_ID ../assembly.scfmap ../assembly.homopolymer-compressed.layout reads.ids\n";
    print CMD "zcat $hifi | seqtk subseq - reads.ids > reads.WORKING.fasta && mv reads.WORKING.fasta reads.fasta\n";
    print CMD "$rsoft/mitohifi/bin/mitohifi -r reads.fasta -f /patched/MitoHiFi/resources/sequence.fasta -g /patched/MitoHiFi/resources/sequence.gb -t \$SLURM_CPUS_PER_TASK -o 2 -a animal -p 90 --rotate-to-ref\n";
    print CMD "minimap2 \\\n";
    print CMD "   -t\$SLURM_CPUS_PER_TASK \\\n";
    print CMD "   --secondary=no \\\n";
    print CMD "   -ax asm5 \\\n";
    print CMD "   $rsoft/mitohifi/resources/sequence.fasta \\\n";
    print CMD "   ./final_mitogenome.fasta \\\n";
    print CMD "   | samtools view -Sb \\\n";
    print CMD "   | samtools sort - \\\n";
    print CMD "   > ${samp}_mito_asm_on_reference.bam \\\n";
    print CMD "samtools index ${samp}_mito_asm_on_reference.bam\n";
    print CMD "bcftools mpileup \\\n";
    print CMD "   -f $rsoft/mitohifi/resources/sequence.fasta \\\n";
    print CMD "   ${samp}_mito_asm_on_reference.bam  \\\n";
    print CMD "  | bcftools call \\\n";
    print CMD "        -mv -Ov \\\n";
    print CMD "        --ploidy 1 \\\n";
    print CMD "        -o ${samp}_mito_asm_on_reference.vcf\n";
    print CMD "fi\n";
    print CMD "\n";
    print CMD "rm -f ./mitohifi.jid\n";
    print CMD "\n";
    print CMD "exit 0\n";
  }

  #  Run, if not finished already.

  if    ($finished)                   { print "$samp/$flav/mitohifi - FINISHED\n"; }
  elsif (-e "$diro/mitohifi.jid")      { print "$samp/$flav/mitohifi - RUNNING\n"; }
  elsif (-e "$diro/mitohifi.FAIL")     { print "$samp/$flav/mitohifi - CRASHED\n"; }
  elsif ($unavail)                    { print "$samp/$flav/mitohifi - UNAVAILABLE ($unavail)\n"; }
  elsif (! $ready)                    { print "$samp/$flav/mitohifi - ASSEMBLY-NOT-READY\n"; }
  elsif (! $$opts{"submit"})          { print "$samp/$flav/mitohifi - READY-TO-COMPUTE\n"; }
  else                                { print "$samp/$flav/mitohifi - SUBMITTED\n"; system("sbatch $diro/mitohifi.sh > $diro/mitohifi.jid"); }
}


1;
