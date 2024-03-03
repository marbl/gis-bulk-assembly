
###############################################################################
#
#  This file is part of a pipeline to (help) run Verkko assemblies on
#  numerous HPRC samples.
#
#  This is a 'United States Government Work', and is released in the public
#  domain.
#
##

package hprc::aws;
require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(awsToLocalPath awsToLocalInfo fetchInfo fetchData);

use strict;
use warnings "all";
no  warnings "uninitialized";

use hprc::samples;


sub awsToLocalPath ($$) {
  my $samp = shift @_;
  my $locf = shift @_;

  $locf =~ s!/!--!g;
  $locf =~ s!^s3:----human-pangenomics--submissions--!aws-data/$samp/!;

  return($locf);
}

sub awsToLocalInfo ($$) {
  my $samp = shift @_;
  my $info = shift @_;

  $info =~ s!/!--!g;
  $info =~ s!^s3:----human-pangenomics--submissions--!aws-info/$samp/!;
  $info .= ".s3ls";

  return($info);
}

#
#  Given a sample name and a file, return the size of the file.
#  Caches directory listings locally in 'aws-info'.
#
sub fetchInfo ($$) {
  my $samp = shift @_;
  my $file = shift @_;
  my $info = awsToLocalInfo($samp, $file);

  if (! -d "aws-info/$samp") {
    system("mkdir -p aws-info/$samp");
  }

  if (! -e "$info") {
    print STDERR "  Fetch AWS ls for $file.\n";
    system("aws --no-sign-request s3 ls $file > $info 2> $info.err");
  }

  if (! -e "$info") {
    die "didn't find '$info'\n";
  }

  my $size = 0;

  open(INF, "< $info") or die "Failed to open '$info' for reading: $!\n";
  while (<INF>) {
    s/^\s+//;
    s/\s+$//;

    my @v = split '\s+', $_;

    $size += $v[2];
  }
  close(INF);

  return $size;
}


#
#  Given a sample name and a file type, fetches all the files into
#  
sub fetchData ($$) {
  my $samp  = shift @_;
  my $type  = shift @_;
  my $files = $samples{$samp}{$type};

  #print STDERR "\n";
  #print STDERR "FETCHING '$type' FOR SAMPLE '$samp'.\n";
  #print STDERR "\n";

  if (! -d "aws-data/$samp") {
    system("mkdir -p aws-data/$samp");
  }

  foreach my $f (@$files) {
    my $awsf = $f;                                     #  So we don't accidentally corrupt the file list.
    my $locf = awsToLocalPath($samp, $f);
    #y $awso = $f =~ s!^s3://human-pangenomics/!!/r;   #  Needed for s3api.

    if (-e $locf) {
      printf "Exists %7s/%-9s %s\n", $samp, "$type:", $locf;
      next;
    }

    printf "Fetch  %7s/%-9s %s\n", $samp, "$type:", $locf;

    #y $c = "aws --no-sign-request s3api get-object --bucket human-pangenomics --key '$awso' --range bytes=0-1048576 '$locf' > $locf.err 2>&1";
    my $c = "aws --no-sign-request s3 cp '$awsf' '$locf' > $locf.err 2>&1";
    my $r = system($c);

    if ($r == 0) {
      next;
    }

    rename '$locf',     '$locf.FAILED';
    rename '$locf.err', '$locf.FAILED.err';

    print "    FAILED: $c\n";
    open(FAIL, "< $locf.err");
    $_ = <FAIL>;
    $_ = <FAIL>;
    while (<FAIL>) {
      print "      $_";
      last;
    }
    close(FAIL);

    print "\n";
  }
}

1;
