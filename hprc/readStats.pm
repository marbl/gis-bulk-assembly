
###############################################################################
#
#  This file is part of a pipeline to (help) run Verkko assemblies on
#  numerous HPRC samples.
#
#  This is a 'United States Government Work', and is released in the public
#  domain.
#
##

package hprc::readStats;
require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(summarizeReads);

use strict;
use warnings "all";
no  warnings "uninitialized";

use hprc::samples;
use hprc::assemble;


sub summarizeReads ($$$) {
  my $samp   = shift @_;
  my $types  = shift @_;
  my $opts   = shift @_;

  foreach my $t (sort keys %$types) {
    my @files = getFiles($samp, $t);

    print "Type $t has ", scalar(@files), " files.\n";

    foreach my $file (@files) {
      my $summ = $file;
      my $name;

      $summ =~ s/.(fasta.gz|fastq.gz|sam|bam|cram)$/.summary/;
      $name =  $summ;
      $name =~ s/^.*--//;
      $name =~ s/.summary$//;

      if (-z $summ) {
        unlink "$summ";
        unlink "$summ.err";
      }

      if (-e $summ) {
        print STDERR " - $name - EXISTS.\n";
      }
      else {
        my $cmd = "$root/software/seqrequester/build/bin/seqrequester summarize '$file' > '$summ' 2> '$summ.err'";

        print STDERR " - $name - ";

        if (system($cmd) != 0) {
          print STDERR "FAILED!\n";
          unlink "$summ.summary";
        }
        else {
          print STDERR "SUCCESS!\n";
          unlink "$summ.err";
        }
      }
    }
  }
}

1;
