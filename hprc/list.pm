
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

use hprc::samples;
use hprc::aws;


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
  my $sample = shift @_;
  my $files  = shift @_;
  my $size   = 0;

  foreach my $file (@$files) {
    $size += fetchInfo($sample, $file);
  }

  return $size / 1024.0 / 1024.0 / 1024.0;
}


sub locsiz ($$) {
  my $sample = shift @_;
  my $files  = shift @_;
  my $size   = 0;

  foreach my $file (@$files) {
    my $locf = awsToLocalPath($sample, $file);
    $size += -s awsToLocalPath($sample, $file);
  }

  return $size / 1024.0 / 1024.0 / 1024.0;
}


sub displayFiles ($$$) {
  my $sample = shift @_;
  my $files  = shift @_;
  my $opts   = shift @_;

  return  if (! exists($$opts{"files"}));

  foreach my $file (@$files) {
    my $locf = awsToLocalPath($sample, $file);
    print "                    - $locf\n";
  }
}


sub displaySample ($$$) {
  my $s      = shift @_;
  my $types  = shift @_;
  my $opts   = shift @_;
  my $awstot = 0;
  my $loctot = 0;

  foreach my $t (keys %$types) {
    $awstot += awssiz($s, $samples{$s}{$t});
    $loctot += locsiz($s, $samples{$s}{$t});
  }

  printf("%7s (%s %s %s/%s %s/%s)\n",
         $s,
         pad( 8, $samples{$s}{'sex'}),
         pad( 5, $samples{$s}{'cohort'}),
         pad(-4, $samples{$s}{'superpop'}),  pad(4, $samples{$s}{'subpop'}),
         pad(-4, $samples{$s}{'prod_year'}),        $samples{$s}{'hifi_prod_site'});
  printf("                AWS / local (GB)\n");
  if (exists($$types{'hifi'}))      { printf("  HiFi:      %6.1f / %6.1f\n", awssiz($s, $samples{$s}{'hifi'}),     locsiz($s, $samples{$s}{'hifi'}));       displayFiles($s, $samples{$s}{'hifi'},     $opts); }
  if (exists($$types{'ont'}))       { printf("  ONT:       %6.1f / %6.1f\n", awssiz($s, $samples{$s}{'ont'}),      locsiz($s, $samples{$s}{'ont'}));        displayFiles($s, $samples{$s}{'ont'},      $opts); }
  if (exists($$types{'hic'}))       { printf("  Hi-C:      %6.1f / %6.1f\n", awssiz($s, $samples{$s}{'hic'}),      locsiz($s, $samples{$s}{'hic'}));        displayFiles($s, $samples{$s}{'hic'},      $opts); }
  if (exists($$types{'ilmn'}))      { printf("  Ilmn:      %6.1f / %6.1f\n", awssiz($s, $samples{$s}{'ilmn'}),     locsiz($s, $samples{$s}{'ilmn'}));       displayFiles($s, $samples{$s}{'ilmn'},     $opts); }
  if (exists($$types{'mat-ilmn'}))  { printf("  Ilmn(mat): %6.1f / %6.1f\n", awssiz($s, $samples{$s}{'mat-ilmn'}), locsiz($s, $samples{$s}{'mat-ilmn'}));   displayFiles($s, $samples{$s}{'mat-ilmn'}, $opts); }
  if (exists($$types{'pat-ilmn'}))  { printf("  Ilmn(pat): %6.1f / %6.1f\n", awssiz($s, $samples{$s}{'pat-ilmn'}), locsiz($s, $samples{$s}{'pat-ilmn'}));   displayFiles($s, $samples{$s}{'pat-ilmn'}, $opts); }
  printf("             %6.1f / %6.1f\n", $awstot, $loctot);
  printf("\n");
}

1;
