
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


sub siz ($$) {
  my $sample = shift @_;
  my $files  = shift @_;
  my $size   = 0;

  foreach my $file (@$files) {
    $size += fetchInfo($sample, $file);
  }

  return $size / 1024.0 / 1024.0 / 1024.0;
}


sub displaySample ($) {
  my $s      = shift @_;
  my $awstot = 0;
  my $loctot = 0;

  printf("%7s (%s %s %s/%s %s/%s)\n",
         $s,
         pad( 8, $samples{$s}{'sex'}),
         pad( 5, $samples{$s}{'cohort'}),
         pad(-4, $samples{$s}{'superpop'}),  pad(4, $samples{$s}{'subpop'}),
         pad(-4, $samples{$s}{'prod_year'}),        $samples{$s}{'hifi_prod_site'});
  printf("                AWS / local (GB)\n");
  printf("  HiFi:      %6.1f / %6.1f\n", siz($s, $samples{$s}{'hifi'}),     0);  $awstot += siz($s, $samples{$s}{'hifi'});
  printf("  ONT:       %6.1f / %6.1f\n", siz($s, $samples{$s}{'ont'}),      0);  $awstot += siz($s, $samples{$s}{'ont'});
  printf("  Hi-C:      %6.1f / %6.1f\n", siz($s, $samples{$s}{'hic'}),      0);  $awstot += siz($s, $samples{$s}{'hic'});
  printf("  Ilmn:      %6.1f / %6.1f\n", siz($s, $samples{$s}{'ilmn'}),     0);  $awstot += siz($s, $samples{$s}{'ilmn'});
  printf("  Ilmn(mat): %6.1f / %6.1f\n", siz($s, $samples{$s}{'mat-ilmn'}), 0);  $awstot += siz($s, $samples{$s}{'mat-ilmn'});
  printf("  Ilmn(pat): %6.1f / %6.1f\n", siz($s, $samples{$s}{'pat-ilmn'}), 0);  $awstot += siz($s, $samples{$s}{'pat-ilmn'});
  printf("             %6.1f / %6.1f\n", $awstot, $loctot);
  printf("\n");
}

1;
