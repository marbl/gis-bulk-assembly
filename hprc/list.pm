
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


sub displaySample ($) {
  my $s      = shift @_;
  my $awstot = 0;
  my $loctot = 0;

  $awstot += awssiz($s, $samples{$s}{'hifi'});
  $awstot += awssiz($s, $samples{$s}{'ont'});
  $awstot += awssiz($s, $samples{$s}{'hic'});
  $awstot += awssiz($s, $samples{$s}{'ilmn'});
  $awstot += awssiz($s, $samples{$s}{'mat-ilmn'});
  $awstot += awssiz($s, $samples{$s}{'pat-ilmn'});

  $loctot += locsiz($s, $samples{$s}{'hifi'});
  $loctot += locsiz($s, $samples{$s}{'ont'});
  $loctot += locsiz($s, $samples{$s}{'hic'});
  $loctot += locsiz($s, $samples{$s}{'ilmn'});
  $loctot += locsiz($s, $samples{$s}{'mat-ilmn'});
  $loctot += locsiz($s, $samples{$s}{'pat-ilmn'});

  printf("%7s (%s %s %s/%s %s/%s)\n",
         $s,
         pad( 8, $samples{$s}{'sex'}),
         pad( 5, $samples{$s}{'cohort'}),
         pad(-4, $samples{$s}{'superpop'}),  pad(4, $samples{$s}{'subpop'}),
         pad(-4, $samples{$s}{'prod_year'}),        $samples{$s}{'hifi_prod_site'});
  printf("                AWS / local (GB)\n");
  printf("  HiFi:      %6.1f / %6.1f\n", awssiz($s, $samples{$s}{'hifi'}),     locsiz($s, $samples{$s}{'hifi'}));
  printf("  ONT:       %6.1f / %6.1f\n", awssiz($s, $samples{$s}{'ont'}),      locsiz($s, $samples{$s}{'ont'}));
  printf("  Hi-C:      %6.1f / %6.1f\n", awssiz($s, $samples{$s}{'hic'}),      locsiz($s, $samples{$s}{'hic'}));
  printf("  Ilmn:      %6.1f / %6.1f\n", awssiz($s, $samples{$s}{'ilmn'}),     locsiz($s, $samples{$s}{'ilmn'}));
  printf("  Ilmn(mat): %6.1f / %6.1f\n", awssiz($s, $samples{$s}{'mat-ilmn'}), locsiz($s, $samples{$s}{'mat-ilmn'}));
  printf("  Ilmn(pat): %6.1f / %6.1f\n", awssiz($s, $samples{$s}{'pat-ilmn'}), locsiz($s, $samples{$s}{'pat-ilmn'}));
  printf("             %6.1f / %6.1f\n", $awstot, $loctot);
  printf("\n");

}

1;
