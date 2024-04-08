
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

  #  Compute?

  foreach my $t (sort keys %$types) {
    my @files = getFiles($samp, $t);
    my $needCompute = 0;

    foreach my $file (@files) {
      my $outs = $file;    $outs =~ s/.(fasta.gz|fastq.gz|sam|bam|cram)$/.summary/;
      my $outl = $file;    $outl =~ s/.(fasta.gz|fastq.gz|sam|bam|cram)$/.readlen/;

      if (-z $outs)   { unlink "$outs"; unlink "$outs.err"; }
      if (-z $outl)   { unlink "$outl"; unlink "$outl.err"; }

      $needCompute = 1   if (! -e $outs);
      $needCompute = 1   if (! -e $outl) && (($t eq "hifi") || ($t eq "ont"));
    }

    next  if ($needCompute == 0);

    print "Computing summaries for type $t with ", scalar(@files), " files.\n";

    foreach my $file (@files) {
      my $outs = $file;    $outs =~ s/.(fasta.gz|fastq.gz|sam|bam|cram)$/.summary/;
      my $outl = $file;    $outl =~ s/.(fasta.gz|fastq.gz|sam|bam|cram)$/.readlen/;
      my $name = $file;    $name =~ s/^.*--//;     #  Display name.
                           $name =~ s/.(fasta.gz|fastq.gz|sam|bam|cram)$//;

      my $cmds = "$root/software/seqrequester/build/bin/seqrequester summarize         '$file' > '$outs.tmp' 2> '$outs.err' && mv '$outs.tmp' '$outs'";
      my $cmdl = "$root/software/seqrequester/build/bin/seqrequester summarize -seqlen '$file' > '$outl.tmp' 2> '$outl.err' && mv '$outl.tmp' '$outl'";

      print " - $name";
      print " -- summary:";
      if (! -e $outs) {
        if (system($cmds) != 0)   { print "FAILED";  unlink "$outs.tmp";  }
        else                      { print "done  ";  unlink "$outs.err";  }
      } else                      { print "exists";                       }

      print " -- readlen:";
      if (($t eq "hic") ||        #  Read length not interesting
          ($t eq "ilmn") ||       #  for illumina data.
          ($t eq "mat-ilmn") ||
          ($t eq "pat-ilmn"))     { print "skipped";                        }
      elsif (! -e $outl) {
        if (system($cmdl) != 0)   { print "FAILED";    unlink "$outl.tmp";  }
        else                      { print "done";      unlink "$outl.err";  }
      } else                      { print "exists";                         }

      print "\n";
    }
  }

  #  Display summaries.

  foreach my $t (sort keys %$types) {
    my @files = getFiles($samp, $t);

    #print "Type $t has ", scalar(@files), " files.\n";

    if ($t eq "hifi") {
      printf "\n";
      printf "               N20 |                N50 |                N80 |               N100 | $samp\n";
      printf " min-len num-reads |  min-len num-reads |  min-len num-reads |  min-len num-reads |\n";
      printf "-------- --------- | -------- --------- | -------- --------- | -------- --------- |\n";
    }

    if ($t eq "ont") {
      printf "\n";
      printf "    0 Kbp < 50 Kbp  |    50 Kbp < 100 Kbp |   100 Kbp < 200 Kbp | 200 Kbp < {max}     | $samp\n";
      printf "     Gbp  num-reads |       Gbp num-reads |       Gbp num-reads |       Gbp num-reads |\n";
      printf "--------- --------- | --------- --------- | --------- --------- | --------- --------- |\n";
    }

    foreach my $file (@files) {
      my $outs = $file;    $outs =~ s/.(fasta.gz|fastq.gz|sam|bam|cram)$/.summary/;
      my $outl = $file;    $outl =~ s/.(fasta.gz|fastq.gz|sam|bam|cram)$/.readlen/;
      my $name = $file;    $name =~ s/^.*--//;     #  Display name.
                           $name =~ s/.(fasta.gz|fastq.gz|sam|bam|cram)$//;

      if ($t eq "hifi") {
        my $len020 = 0;  my $idx020 = 0;
        my $len050 = 0;  my $idx050 = 0;
        my $len080 = 0;  my $idx080 = 0;
        my $len100 = 0;  my $idx100 = 0;

        open(SUM, "< $outs");
        while (<SUM>) {
          if (m/^00020\s+(\d+)\s+(\d+)\s+/)   { $len020 = $1;  $idx020 = $2; }
          if (m/^00050\s+(\d+)\s+(\d+)\s+/)   { $len050 = $1;  $idx050 = $2; }
          if (m/^00080\s+(\d+)\s+(\d+)\s+/)   { $len080 = $1;  $idx080 = $2; }
          if (m/^00100\s+(\d+)\s+(\d+)\s+/)   { $len100 = $1;  $idx100 = $2; }
        }

        close(SUM);

        printf "%8d %9d | %8d %9d | %8d %9d | %8d %9d | %s\n",
            $len020, $idx020, $len050, $idx050, $len080, $idx080, $len100, $idx100, $name;
      }

      if ($t eq "ont") {
        my $n000 = 0;  my $s000 = 0;
        my $n050 = 0;  my $s050 = 0;
        my $n100 = 0;  my $s100 = 0;
        my $n200 = 0;  my $s200 = 0;
        my $lmax = 0;

        open(LEN, "< $outl");
        while (<LEN>) {
          if (m/^\d+\s+(\d+)\s+/) {
            if    ($1 >= 200000)  { $n200++;  $s200 += $1; }
            elsif ($1 >= 100000)  { $n100++;  $s100 += $1; }
            elsif ($1 >=  50000)  { $n050++;  $s050 += $1; }
            else                  { $n000++;  $s000 += $1; }

            $lmax = $1   if ($lmax < $1);
          }
        }
        close(LEN);

        printf "%9.3f %9d | %9.3f %9d | %9.3f %9d | %9.3f %9d | %s\n",
            $s000 / 1024.0 / 1024.0 / 1024.0, $n000,
            $s050 / 1024.0 / 1024.0 / 1024.0, $n050,
            $s100 / 1024.0 / 1024.0 / 1024.0, $n100,
            $s200 / 1024.0 / 1024.0 / 1024.0, $n200,
            $name;
      }
    }
  }
}

1;
