
###############################################################################
#
#  This file is part of a pipeline to (help) run Verkko assemblies on
#  numerous HPRC samples.
#
#  This is a 'United States Government Work', and is released in the public
#  domain.
#
##

package hprc::list;
require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(displaySample);

use strict;
use warnings "all";
no  warnings "uninitialized";
use experimental "for_list";

use hprc::samples;
use hprc::aws;



my $displaySummaryIndex = 0;

sub displaySummary ($$$) {
  my $samp   = shift @_;
  my $types  = shift @_;
  my $opts   = shift @_;

  my ($HnReads, $OnReads, $TnReads, $CnReads) = ( 0, 0, 0, 0 );
  my ($HnBases, $OnBases, $TnBases, $CnBases) = ( 0, 0, 0, 0 );
  my ($isHic12, $isHic) = (0, 0);

  if ($displaySummaryIndex == 0) {
    print "--------  ------- HiFi -------  -------- ONT -------  --- Child --  -------------- Trio --------------  -------- Hi-C -------\n";
    print "                                                                                  maternal   paternal\n";
    print "sample    coverage   het'osity  coverage   het'osity     het'osity  coverage   het'osity    het'osity  coverage   het'osity\n";
    print "--------  --------------------  --------------------  ------------  ---------------------------------  --------------------\n";
    $displaySummaryIndex = 1;
  }

  #  Read sequence summaries.

  foreach my $t (sort keys %$types) {
    my %files = getFileMap($samp, $t, "even-those-that-don't-exist");
    my $nReads = 0;
    my $nBases = 0;

    foreach my $file (values %files) {
      my $outs = $file;    $outs =~ s/.(fasta.gz|fastq.gz|sam|bam|cram)$/.summary/;
      my $outl = $file;    $outl =~ s/.(fasta.gz|fastq.gz|sam|bam|cram)$/.readlen/;
      my $name = $file;    $name =~ s/^.*--//;     #  Display name.
      $name =~ s/.(fasta.gz|fastq.gz|sam|bam|cram)$//;

      if (-e "$outs") {
        open(SUM, "< $outs") or die "Failed to open '$outs': $!\n";
        while (<SUM>) {
          if (m/^00100\s+(\d+)\s+(\d+)\s+(\d+)\s+/) {
            $nReads += $2;
            $nBases += $3;
          }
        }
        close(SUM);
      }
    }

    if ($t eq "hifi")          { $HnReads += $nReads;  $HnBases += $nBases; }
    if ($t eq "hifi-cutadapt") { $HnReads += $nReads;  $HnBases += $nBases; }
    if ($t eq "ont")           { $OnReads += $nReads;  $OnBases += $nBases; }
    if ($t eq "mat-ilmn")      { $TnReads += $nReads;  $TnBases += $nBases; }
    if ($t eq "pat-ilmn")      { $TnReads += $nReads;  $TnBases += $nBases; }
    if ($t eq "hic1" || $t eq "hic2") {
       $CnReads += $nReads;  $CnBases += $nBases;
       print STDERR "Warning: have both hic1/2 and hic files so coverage will be inaccurate\n" if ($isHic);
       $isHic12=1;
    }
    if ($t eq "hic")  {
       $CnReads += $nReads;  $CnBases += $nBases;
       print STDERR "Warning: have both hic1/2 and hic files so coverage will be inaccurate\n" if ($isHic12);
       $isHic=1;
    }
  }

  #  Read genomescope summary.

  my ($HhetMin, $HhetMax) = ( 0, 1 );
  my ($ChetMin, $ChetMax) = ( 0, 1 );
  my ($MhetMin, $MhetMax) = ( 0, 1 );
  my ($PhetMin, $PhetMax) = ( 0, 1 );

  if (-e "$root/hprc-data/$samp/genomescope/hifi/summary.txt") {
    open(GS, "< $root/hprc-data/$samp/genomescope/hifi/summary.txt") or die "Failed to open 'hprc-data/$samp/genomescope/hifi/summary.txt' for reading: $!\n";
    while (<GS>) {
      if (m/Heterozygous\s+\(ab\)\s+(\d+.\d+)%\s+(\d+.\d+)%/) {
        $HhetMin = $1;
        $HhetMax = $2;
      }
    }
    close(GS);
  }
  if (-e "$root/hprc-data/$samp/genomescope/ilmn/summary.txt") {
    open(GS, "< $root/hprc-data/$samp/genomescope/ilmn/summary.txt") or die "Failed to open 'hprc-data/$samp/genomescope/ilmn/summary.txt' for reading: $!\n";
    while (<GS>) {
      if (m/Heterozygous\s+\(ab\)\s+(\d+.\d+)%\s+(\d+.\d+)%/) {
        $ChetMin = $1;
        $ChetMax = $2;
      }
    }
    close(GS);
  }
  if (-e "$root/hprc-data/$samp/genomescope/mati/summary.txt") {
    open(GS, "< $root/hprc-data/$samp/genomescope/mati/summary.txt") or die "Failed to open 'hprc-data/$samp/genomescope/mati/summary.txt' for reading: $!\n";
    while (<GS>) {
      if (m/Heterozygous\s+\(ab\)\s+(\d+.\d+)%\s+(\d+.\d+)%/) {
        $MhetMin = $1;
        $MhetMax = $2;
      }
    }
    close(GS);
  }
  if (-e "$root/hprc-data/$samp/genomescope/pati/summary.txt") {
    open(GS, "< $root/hprc-data/$samp/genomescope/pati/summary.txt") or die "Failed to open 'hprc-data/$samp/genomescope/pati/summary.txt' for reading: $!\n";
    while (<GS>) {
      if (m/Heterozygous\s+\(ab\)\s+(\d+.\d+)%\s+(\d+.\d+)%/) {
        $PhetMin = $1;
        $PhetMax = $2;
      }
    }
    close(GS);
  }

  sub showHet ($$) {
    my $hmin = shift @_;
    my $hmax = shift @_;
    my $hets = "";

    if (($hmin <= 0.0) || ($hmax >= 1.0)) {
      $hets = "            ";
    }
    else {
      $hets = sprintf("%5.3f-%5.3f%%", $hmin, $hmax);
    }
  }

  printf "%-8s  %6.2fx %12s  %6.2fx %12s  %12s  %6.2fx %12s %12s  %6.2fx %12s\n",
      $samp,
      $HnBases / 3100000000.0, showHet($HhetMin, $HhetMax),                                #  HiFi
      $OnBases / 3100000000.0, showHet(   0.000,    1.000),                                #  ONT
                               showHet($ChetMin, $ChetMax),                                #  Child
      $TnBases / 3100000000.0, showHet($MhetMin, $MhetMax), showHet($PhetMin, $PhetMax),   #  Parents
      $CnBases / 3100000000.0, showHet(   0.000,    1.000);                                #  Hi-C
}


sub displayDetails ($$$) {
  my $samp  = shift @_;
  my $type  = shift @_;
  my $opts  = shift @_;
  my %files = getFileMap($samp, $type, "even-those-that-don't-exist");

  #print "Type $t has ", scalar(values %files), " files.\n";

  if ($type eq "hifi" || $type eq "hifi-cutadapt") {
    printf "\n";
    printf "               N20 |                N50 |                N80 |               N100 | $samp\n";
    printf " min-len num-reads |  min-len num-reads |  min-len num-reads |  min-len num-reads |\n";
    printf "-------- --------- | -------- --------- | -------- --------- | -------- --------- |\n";
  }

  if ($type eq "ont") {
    printf "\n";
    printf "    0 Kbp < 50 Kbp  |    50 Kbp < 100 Kbp |   100 Kbp < 200 Kbp | 200 Kbp < {max}     | $samp\n";
    printf "     Gbp  num-reads |       Gbp num-reads |       Gbp num-reads |       Gbp num-reads |\n";
    printf "--------- --------- | --------- --------- | --------- --------- | --------- --------- |\n";
  }

  my ($hifiR020, $hifiB020) = (0, 0);
  my ($hifiR050, $hifiB050) = (0, 0);
  my ($hifiR080, $hifiB080) = (0, 0);
  my ($hifiR100, $hifiB100) = (0, 0);

  my ($ontR200,  $ontB200)  = (0, 0);
  my ($ontR100,  $ontB100)  = (0, 0);
  my ($ontR050,  $ontB050)  = (0, 0);
  my ($ontR000,  $ontB000)  = (0, 0);

  foreach my $file (values %files) {
    my $outs = $file;    $outs =~ s/.(fasta.gz|fastq.gz|sam|bam|cram)$/.summary/;
    my $outl = $file;    $outl =~ s/.(fasta.gz|fastq.gz|sam|bam|cram)$/.readlen/;
    my $name = $file;    $name =~ s/^.*--//;     #  Display name.
    $name =~ s/.(fasta.gz|fastq.gz|sam|bam|cram)$//;

    if ($type eq "hifi" || $type eq "hifi-cutadapt") {
      my ($len020, $idx020) = (0, 0);
      my ($len050, $idx050) = (0, 0);
      my ($len080, $idx080) = (0, 0);
      my ($len100, $idx100) = (0, 0);

      open(SUM, "< $outs");
      while (<SUM>) {
        if (m/^00020\s+(\d+)\s+(\d+)\s+(\d+)\s+/)   { $len020 = $1;  $idx020 = $2;  $hifiR020 += $2;  $hifiB020 += $3; }
        if (m/^00050\s+(\d+)\s+(\d+)\s+(\d+)\s+/)   { $len050 = $1;  $idx050 = $2;  $hifiR050 += $2;  $hifiB050 += $3; }
        if (m/^00080\s+(\d+)\s+(\d+)\s+(\d+)\s+/)   { $len080 = $1;  $idx080 = $2;  $hifiR080 += $2;  $hifiB080 += $3; }
        if (m/^00100\s+(\d+)\s+(\d+)\s+(\d+)\s+/)   { $len100 = $1;  $idx100 = $2;  $hifiR100 += $2;  $hifiB100 += $3; }
      }

      close(SUM);

      printf "%8d %9d | %8d %9d | %8d %9d | %8d %9d | %s\n",
          $len020, $idx020, $len050, $idx050, $len080, $idx080, $len100, $idx100, $name;
    }

    if ($type eq "ont") {
      my ($n000, $s000) = (0, 0);
      my ($n050, $s050) = (0, 0);
      my ($n100, $s100) = (0, 0);
      my ($n200, $s200) = (0, 0);
      my  $lmax = 0;

      open(LEN, "< $outl");
      while (<LEN>) {
        if (m/^\d+\s+(\d+)\s+/) {
          if    ($1 >= 200000)  { $n200++;  $s200 += $1; }
          elsif ($1 >= 100000)  { $n100++;  $s100 += $1; }
          elsif ($1 >=  50000)  { $n050++;  $s050 += $1; }
          else                  { $n000++;  $s000 += $1; }

          if    ($1 >= 200000)  { $ontR200++;  $ontB200 += $1; }
          if    ($1 >= 100000)  { $ontR100++;  $ontB100 += $1; }
          if    ($1 >=  50000)  { $ontR050++;  $ontB050 += $1; }
          { $ontR000++;  $ontB000 += $1; }

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

  #  Print type summaries

  if ($type eq "hifi" || $type eq "hifi-cutadapt") {
    printf "-------- --------- | -------- --------- | -------- --------- | -------- --------- |\n";
    printf "%7.3fx %9d | %7.3fx %9d | %7.3fx %9d | %7.3fx %9d |\n",
        $hifiB020 / 3100000000, $hifiR020, $hifiB050 / 3100000000, $hifiR050, $hifiB080 / 3100000000, $hifiR080, $hifiB100 / 3100000000, $hifiR100;
  }

  if ($type eq "ont") {
    printf "--------- --------- | --------- --------- | --------- --------- | --------- --------- |\n";
    printf "%8.3fx %9d | %8.3fx %9d | %8.3fx %9d | %8.3fx %9d |\n",
        $ontB000 / 3100000000, $ontR000, $ontB050 / 3100000000, $ontR050, $ontB100 / 3100000000, $ontR100, $ontB200 / 3100000000, $ontR200;
  }
}




sub pad ($$) {
  my $n = shift @_;
  my $s = shift @_;
  my $S = " " x abs($n);
  if ($n < 0) {
    #print "NEG: '$S$s' $n - '", substr("$S$s", $n), "'\n";
    return substr("$S$s", $n);
  } else {
    #print "POS: '$s;$S' $n - '", substr("$s;$S", 0, $n), "'\n";
    return substr("$s;$S", 0, $n);
  }
}


sub awssiz ($$) {
  my $samp  = shift @_;
  my $files = shift @_;
  my $size  = 0;

  foreach my $file (keys %$files) {
    $size += getRemoteSize($samp, $file);
  }

  return $size;
}


sub locsiz ($$) {
  my $samp  = shift @_;
  my $files = shift @_;
  my $size  = 0;

  foreach my $file (values %$files) {
    $size += getLocalSize($samp, $file);
  }

  return $size;
}


sub displayFiles ($$$$) {
  my $samp    = shift @_;
  my $type    = shift @_;
  my $typen   = shift @_;
  my $opts    = shift @_;
  my %filemap = getFileMap($samp, $type, "and-the-one-that-have-been-deleted");

  my $awss = awssiz($samp, \%filemap);
  my $locs = locsiz($samp, \%filemap);

  $awss = 0  if ($type eq "hifi-cutadapt");

  if (exists($$opts{"files"})) {
    foreach my ($awsf, $locf) (%filemap) {
      my $r = getRemoteSize($samp, $awsf);
      my $l = getLocalSize($samp, $locf);
      my $z = ($l != $r) ? "!" : " ";

      printf("%7s/%-15s  %6.1f %s %6.1f  %s\n", $samp, $type, $r, $z, $l, $locf);
    }
    printf("%7s/%-15s  %6.1f   %6.1f\n", $samp, $type, $awss, $locs);
  }
  else {
    printf("%7s/%-15s  %6.1f   %6.1f\n", $samp, $type, $awss, $locs);
  }
}


sub displaySizes ($$$) {
  my $samp   = shift @_;
  my $types  = shift @_;
  my $opts   = shift @_;
  my $awstot = 0;
  my $loctot = 0;

  #  If 'hifi' is requested, also show 'hifi-cutadapt'.
  $$types{'hifi-cutadapt'} = 1   if ($$types{'hifi'});

  foreach my $t (keys %$types) {         #  Sum aws and loc sizes for any types supplied, but if 'hifi-cutadapt'
    my %fm = getFileMap($samp, $t, 1);   #  is explictly requested, ignore the aws size (it's included in 'hifi').

    $awstot += awssiz($samp, \%fm)   if ($t ne "hifi-cutadapt");   #  See displayFiles() above too!
    $loctot += locsiz($samp, \%fm);
  }

  printf("%7s (%s %s %s/%s %s/%s)\n",
         $samp,
         pad( 8, $samples{$samp}{'sex'}),
         pad( 5, $samples{$samp}{'cohort'}),
         pad(-4, $samples{$samp}{'superpop'}),  pad(4, $samples{$samp}{'subpop'}),
         pad(-4, $samples{$samp}{'prod_year'}),        $samples{$samp}{'hifi_prod_site'});

  if (exists($$opts{"files"})) {
    print "-----------------------  -- SIZE (GB) --  ----------------------------------------\n";
    print "sample  read-type           aws    local  path\n";
    print "------- ---------------  ------ - ------  ----------\n";
  }
  else {
    print "-----------------------  -- SIZE (GB) --\n";
    print "sample  read-type           aws    local\n";
    print "------- ---------------  ------ - ------\n";
  }

  if (exists($$types{'hifi'}))          { displayFiles($samp, 'hifi',          "HiFi:           ", $opts); }
  if (exists($$types{'hifi-cutadapt'})) { displayFiles($samp, 'hifi-cutadapt', "HiFi(cutadapt): ", $opts); }
  if (exists($$types{'ont'}))           { displayFiles($samp, 'ont',           "ONT:            ", $opts); }
  if (exists($$types{'hic'}))           { displayFiles($samp, 'hic',           "Hi-C:           ", $opts); }
  if (exists($$types{'ilmn'}))          { displayFiles($samp, 'ilmn',          "Ilmn:           ", $opts); }
  if (exists($$types{'mat-ilmn'}))      { displayFiles($samp, 'mat-ilmn',      "Ilmn(mat):      ", $opts); }
  if (exists($$types{'pat-ilmn'}))      { displayFiles($samp, 'pat-ilmn',      "Ilmn(pat):      ", $opts); }

  printf("%7s                  %6.1f   %6.1f\n", $samp, $awstot, $loctot);
  printf("\n");
}



sub displaySample ($$$) {
  my $samp   = shift @_;
  my $types  = shift @_;
  my $opts   = shift @_;

  if    ($$opts{'summary'} == 1) {
    displaySummary($samp, $types, $opts);
  }

  elsif ($$opts{'details'} == 1) {
    foreach my $type (sort keys %$types) {
      displayDetails($samp, $type, $opts);
    }
  }

  else {
    displaySizes($samp, $types, $opts);
  }
}

1;
