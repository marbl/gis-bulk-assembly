#!/usr/bin/perl

use strict;
use List::Util qw(min max);

my %ctgcov;
my %refcov;
my %minp;
my %maxp;
my %len;
my %assign;

#
#  Expects mashmap of assembly-vs-chromosome on stdin.
#
#      mashmap \
#        --ref /data/Phillippy/references/T2T-CHM13/chm13v2.0.fa \
#        --query ../assembly.fasta \
#        --perc_identity 95 \
#        --segLength 100000 \
#        --threads 4 \
#        --output assembly-chromosome.mashmap.out
#

while (<>) {
  #my @v = split '\s', $_;

  my ($qn, $qlen, $qbgn, $qend, $ori, $rn, $rlen, $rbgn, $rend, $rm, $bl, $mq, $id, $kc) = split '\s', $_;

  if ($id =~ m/^id:f:/)   {  $id =~ s/id:f://;  $id *= 100.0;  }   #  Estimated ANI
  if ($kc =~ m/^kc:f:/)   {  $kc =~ s/kc:f://;                 }   #  Kmer complexity, >= v3.0.6

  $ori = "fwd"  if ($ori eq "+");
  $ori = "rev"  if ($ori eq "-");

  my $tag    = "$qn\t$rn\t$ori";

  next if ($qend - $qbgn < 500000);
  next if ($id           < 99.0);

  $ctgcov{$tag} += ($qend - $qbgn) / $qlen;
  $refcov{$tag} += ($rend - $rbgn) / $rlen;

  $minp{$tag} = 999999999   if (!exists($minp{$tag}));
  $maxp{$tag} = 0           if (!exists($maxp{$tag}));

  $minp{$tag} = min($rbgn, $minp{$tag});
  $maxp{$tag} = max($rend, $maxp{$tag});

  $len{$qn} = $qlen;
}

foreach my $k (sort keys %ctgcov) {
  my ($ctg, $chr, $ori) = split '\t', $k;

  if ($ctgcov{$k} > 0.1) {   #  Count if a contig is assigned to more than
    $assign{$ctg}++;         #  one chr/ori.
  }
}


print "-------ASSEMBLY-CONTIG--------   --------REFERENCE-CHROMOSOME--------  ------- -------\n";
print "contig-name         length ori   chromosome-name      begin end        ctg-cov ref-cov\n";
print "---------------- --------- ---   ---------------- --------- ---------  ------- -------\n";

sub byChrPos {
  my ($actg, $achr, $aori, $achri) = split '\t', $a;
  my ($bctg, $bchr, $bori, $bchri) = split '\t', $b;

  my $actgp = $actg;   $actgp = $1  if ($actgp =~ m/^(.*?)\d+$/);
  my $bctgp = $bctg;   $bctgp = $1  if ($bctgp =~ m/^(.*?)\d+$/);

  return -1  if ($actgp lt $bctgp);       #  Sort by contig class first.
  return  1  if ($actgp gt $bctgp);

  $achri = $1   if ($achr =~ m/(\d+)/);   #  Then by chromosome number or name.
  $bchri = $1   if ($bchr =~ m/(\d+)/);

  if (defined($achri) && defined($bchri)) {
    return -1  if ($achri < $bchri);
    return  1  if ($achri > $bchri);
  }
  else {
    return -1  if ($achr lt $bchr);
    return  1  if ($achr gt $bchr);
  }

  my $abgn = $minp{$a};                   #  Then by chromosome position.
  my $bbgn = $minp{$b};

  return $abgn <=> $bbgn;
}


foreach my $k (sort byChrPos (keys %ctgcov)) {
  my ($ctg, $chr, $ori) = split '\t', $k;

  if ($ctgcov{$k} > 0.1) {
    my $min = $minp{$k};
    my $max = $maxp{$k};
    my $len = $len{$ctg};

    printf "%-16s %9d %3s  %-16s %9d-%-9d %8.3f %8.3f%s\n",
        $ctg, $len, $ori, $chr, $min, $max, $ctgcov{$k}, $refcov{$k}, ($assign{$ctg} == 1) ? "" : " SPLIT";
  }
}

exit(0);
