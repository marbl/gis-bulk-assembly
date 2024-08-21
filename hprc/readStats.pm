
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
use hprc::aws;

sub summarizeReads ($$$) {
  my $samp   = shift @_;
  my $types  = shift @_;
  my $opts   = shift @_;

  foreach my $type (sort keys %$types) {
    my %files = getFileMap($samp, $type, "include-files-that-don't-exist");
    my $needCompute = 0;

    foreach my $file (values %files) {
      my $outs = $file;    $outs =~ s/.(fasta.gz|fastq.gz|sam|bam|cram)$/.summary/;
      my $outl = $file;    $outl =~ s/.(fasta.gz|fastq.gz|sam|bam|cram)$/.readlen/;
      my $name = $file;    $name =~ s/.(fasta.gz|fastq.gz|sam|bam|cram)$//;  $name =~ s!/!--!;  $name =~ s!^.*--!!;

      if (-z $outs)   { unlink "$outs"; unlink "$outs.err"; }
      if (-z $outl)   { unlink "$outl"; unlink "$outl.err"; }

      my $needs = (! -e $outs);
      my $needl = (! -e $outl) && (($type eq "hifi") || ($type eq "ont") || ($type eq "hifi-cutadapt"));

      printf "%7s/%-9s - summary:%-7s readlen:%-7s - %s\n", $samp, $type, $needs ? "missing" : "ok", $needl ? "missing" : "ok", $name;

      $needCompute = $needs || $needl;
    }

    if ($needCompute) {
      print STDERR "Computing summaries for type $type with ", scalar(values %files), " files.\n";

      foreach my $aws (keys %files) {
        print STDERR "$aws -> $files{$aws}\n";
      }

      foreach my $file (values %files) {
        my $outs = $file;    $outs =~ s/.(fasta.gz|fastq.gz|sam|bam|cram)$/.summary/;
        my $outl = $file;    $outl =~ s/.(fasta.gz|fastq.gz|sam|bam|cram)$/.readlen/;
        my $name = $file;    $name =~ s/.(fasta.gz|fastq.gz|sam|bam|cram)$//;  $name =~ s!/!--!;  $name =~ s!^.*--!!;

        my $cmds = "$root/software/seqrequester/build/bin/seqrequester summarize         '$file' > '$outs.tmp' 2> '$outs.err' && mv '$outs.tmp' '$outs'";
        my $cmdl = "$root/software/seqrequester/build/bin/seqrequester summarize -seqlen '$file' > '$outl.tmp' 2> '$outl.err' && mv '$outl.tmp' '$outl'";

        print STDERR "   s $cmds\n";
        print STDERR "   l $cmdl\n";
        print STDERR " - $name";
        print STDERR " -- summary:";
        if (! -e $outs) {
          if    (! -e $file)          { print STDERR  "MISSING";                       }
          elsif (system($cmds) != 0)  { print STDERR  "FAILED ";  unlink "$outs.tmp";  }
          else                        { print STDERR  "done   ";  unlink "$outs.err";  }
        } else                        { print STDERR  "exists ";                       }

        print STDERR  " -- readlen:";
        if (($type eq "hic") ||        #  Read length not interesting
            ($type eq "ilmn") ||       #  for illumina data.
            ($type eq "mat-ilmn") ||
            ($type eq "pat-ilmn"))     { print STDERR  "skipped";                      }
        elsif (! -e $outl) {
          if    (! -e $file)          { print STDERR  "MISSING";                       }
          elsif (system($cmdl) != 0)  { print STDERR  "FAILED";   unlink "$outl.tmp";  }
          else                        { print STDERR  "done";     unlink "$outl.err";  }
        } else                        { print STDERR  "exists";                        }

        print STDERR  "\n";
      }
    }
  }
}

1;
