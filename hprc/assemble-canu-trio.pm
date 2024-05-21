
###############################################################################
#
#  This file is part of a pipeline to (help) run Verkko assemblies on
#  numerous HPRC samples.
#
#  This is a 'United States Government Work', and is released in the public
#  domain.
#
##


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

1;
