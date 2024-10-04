
###############################################################################
#
#  This file is part of a pipeline to (help) run Verkko assemblies on
#  numerous HPRC samples.
#
#  This is a 'United States Government Work', and is released in the public
#  domain.
#
##

package hprc::assemble;
require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(isFinished cleanupAssembly archiveAssembly computeAssembly);

use strict;
use warnings "all";
no  warnings "uninitialized";

use hprc::samples;
use hprc::aws;

require "hprc/assemble-status.pm";

require "hprc/assemble-canu-hifi.pm";
require "hprc/assemble-canu-trio.pm";

require "hprc/assemble-verkko-base.pm";
require "hprc/assemble-verkko-trio.pm";
require "hprc/assemble-verkko-hi-c.pm";
require "hprc/assemble-verkko-thic.pm";

require "hprc/assemble-cleanup.pm";
require "hprc/assemble-archive.pm";

sub getScreenOption($) {
   my $verkko = shift @_;
   my $isNewScreen = `$verkko/bin/verkko --help |grep -c screen-human`;

   if ($isNewScreen != 0) {
      return "--screen-human-contaminants";
   } else {
      return "--screen human";
   }
}

sub submitIf ($$$$$$$$) {
  my $samp    = shift @_;
  my $flav    = shift @_;
  my $scr     = shift @_;
  my $compl   = shift @_;
  my $previd  = shift @_;
  my $missing = shift @_;
  my $unavail = shift @_;
  my $submit  = shift @_;
  my $jid     = undef;
  my $wait    = ($previd == 0) ? "" : "--depend=afterany:$previd";

  if    ($compl)              { print "$samp/$flav - FINISHED\n";                                return undef; }
  elsif (-e "$scr.jid")       { print "$samp/$flav - RUNNING\n";                                               }
  elsif (-e "$scr.err")       { print "$samp/$flav - CRASHED\n";                                 return undef; }
  elsif ($unavail)            { print "$samp/$flav - UNAVAILABLE ($unavail)\n";                  return undef; }
  elsif ($missing)            { print "$samp/$flav - MISSING-INPUTS (missing $missing)\n";       return undef; }
  elsif (! $submit)           { print "$samp/$flav - READY-TO-COMPUTE\n";                        return undef; }
  else                        { print "$samp/$flav - SUBMITTED\n";  system("sbatch $wait $scr.sh > $scr.jid"); }

  #  An existing scr.jid means we either just submitted, or the assembly is
  #  already running.  Return that jobid for later use.

  if (-e "$scr.jid") {
    open(JID, "< $scr.jid") or die "Failed to open '$scr.jid' for reading: $!\n";
    $jid = int(<JID>);
    close(JID);
  }

  return $jid;
}


sub computeAssembly ($$) {
  my $samp    = shift @_;
  my $opts    = shift @_;
  my $flav    = $$opts{"flavor"};
  my $submit  = $$opts{"submit"};

  #  Check which, if any, outputs exists.

  my $complB = isFinished($samp, "verkko-base");
  my $complT = isFinished($samp, "verkko-trio");
  my $complH = isFinished($samp, "verkko-hi-c");
  my $complC = isFinished($samp, "verkko-thic");

  my $compl  = isFinished($samp, $flav);

  #  Fail aggressively if the base assembly isn't present, but others are.

  print STDERR "$samp/$flav - ERROR: a Trio assembly exists but no Base assembly does!  Fix manually.\n"   if (!$complB) && ($complT);
  print STDERR "$samp/$flav - ERROR: a Hi-C assembly exists but no Base assembly does!  Fix manually.\n"   if (!$complB) && ($complH);
  print STDERR "$samp/$flav - ERROR: a THiC assembly exists but no Base assembly does!  Fix manually.\n"   if (!$complB) && ($complC);
  print STDERR "$samp/$flav - ERROR: a THiC assembly exists but no Trio assembly does!  Fix manually.\n"   if (!$complT) && ($complC);

  return   if (((!$complB) && ($complT)) ||
               ((!$complB) && ($complH)) ||
               ((!$complB) && ($complC)) ||
               ((!$complT) && ($complC)));

  #  Check which inputs exist.

  my $hifi = getDownloadedFiles($samp, "hifi-cutadapt");   #  Return cutadapt form of hifi data.
  my $nano = getDownloadedFiles($samp, "ont");

  my $mati = getDownloadedFiles($samp, "mat-ilmn");
  my $pati = getDownloadedFiles($samp, "pat-ilmn");

  my $hic1 = getDownloadedFiles($samp, "hic1");
  my $hic2 = getDownloadedFiles($samp, "hic2");

  my $hifiMissing   =  (($hifi eq "") && (numFiles($samp, "hifi")     > 0));
  my $nanoMissing   =  (($nano eq "") && (numFiles($samp, "ont")      > 0));
  my $trioMissing   = ((($mati eq "") && (numFiles($samp, "mat-ilmn") > 0)) ||
                       (($pati eq "") && (numFiles($samp, "pat-ilmn") > 0)));
  my $hicMissing    = ((($hic1 eq "") && (numFiles($samp, "hic")      > 0)) ||
                       (($hic2 eq "") && (numFiles($samp, "hic")      > 0)));
  my $hapmerMissing = locateHapmers($samp);

  my $tu  = ((numFiles($samp, "mat-ilmn") == 0) || (numFiles($samp, "pat-ilmn") == 0));
  my $hu  =  (numFiles($samp, "hic")      == 0);

#  my $unavailT = "";
#  my $unavailH = "";
#  my $unavailC = "";
#
#  if    ($tu && $hu) {
#    $unavailT = "no trio data";
#    $unavailH = "no hi-c data";
#    $unavailC = "no trio or hi-c data";
#  }
#  elsif ($tu) {
#    $unavailT = "no trio data";
#    $unavailC = "no trio data";
#  }
#  elsif ($hu) {
#    $unavailH = "no hi-c data";
#    $unavailC = "no hi-c data";
#  }

  my ($unavailT, $unavailH, $unavailC) = dataAvailable($samp);

  my $baseAsmMissing =   (! -e "$rasm/$samp/verkko-base/contigs.fasta");
  my $trioAsmMissing =   (! -e "$rasm/$samp/verkko-trio/assembly.fasta");

  #  Fail if all inputs are not created or downloaded.

  my $missingCH;
  my $missingCT;
  my $missingB;
  my $missingT;
  my $missingH;
  my $missingC;
  my $missing;

  my $iflavor;

  sub combine (@) {
    if (scalar(@_) == 0)   { return undef; }
    elsif (scalar(@_) > 1) { return (join ', ', @_[0..@_-2]) . " and $_[-1]"; }
    else                   { return pop @_; }
  }

  if  ($flav eq "canu-hifi") {
    my @m;
    push @m, "hifi-cutadapt"     if ($hifiMissing);
    $missingCH = combine(@m);
  }
  if  ($flav eq "canu-trio") {
    my @m;
    push @m, "ont"               if ($nanoMissing);
    push @m, "trio"              if ($trioMissing);
    $missingCT = combine(@m);
  }
  if (($flav eq "verkko-base") || ($flav eq "verkko-full")) {
    my @m;
    push @m, "hifi-cutadapt"     if ($hifiMissing);
    push @m, "ont"               if ($nanoMissing);
    #ush @m, "trio"              if ($trioMissing);
    #ush @m, "hic"               if ($hicMissing);
    #ush @m, "base assembly"     if ($baseAsmMissing);
    #ush @m, "trio assembly"     if ($trioAsmMissing);
    #ush @m, "hapmer databases"  if ($hapmerMissing);
    $missingB = combine(@m);
  }
  if (($flav eq "verkko-trio") || ($flav eq "verkko-full")) {
    my @m;
    push @m, "hifi-cutadapt"     if ($hifiMissing);
    push @m, "ont"               if ($nanoMissing);
    push @m, "trio"              if ($trioMissing);
    #ush @m, "hic"               if ($hicMissing);
    push @m, "base assembly"     if ($baseAsmMissing) && ($flav ne "verkko-full");
    #ush @m, "trio assembly"     if ($trioAsmMissing);
    push @m, "hapmer databases"  if ($hapmerMissing);
    $missingT = combine(@m);
  }
  if (($flav eq "verkko-hi-c") || ($flav eq "verkko-full")) {
    my @m;
    push @m, "hifi-cutadapt"     if ($hifiMissing);
    push @m, "ont"               if ($nanoMissing);
    #ush @m, "trio"              if ($trioMissing);
    push @m, "hic"               if ($hicMissing);
    push @m, "base assembly"     if ($baseAsmMissing) && ($flav ne "verkko-full");
    #ush @m, "trio assembly"     if ($trioAsmMissing);
    #ush @m, "hapmer databases"  if ($hapmerMissing);
    $missingH = combine(@m);
  }
  if (($flav eq "verkko-thic") || ($flav eq "verkko-full")) {
    my @m;
    push @m, "hifi-cutadapt"     if ($hifiMissing);
    push @m, "ont"               if ($nanoMissing);
    push @m, "trio"              if ($trioMissing);
    push @m, "hic"               if ($hicMissing);
    push @m, "base assembly"     if ($baseAsmMissing) && ($flav ne "verkko-full");
    push @m, "trio assembly"     if ($trioAsmMissing) && ($flav ne "verkko-full");
    push @m, "hapmer databases"  if ($hapmerMissing);
    $missingC = combine(@m);
  }

  if ($flav eq "") {
    $iflavor  = 1;
    $flav     = "{not-specified}";
  }

  if (($flav ne "canu-hifi")   &&
      ($flav ne "canu-trio")   &&
      ($flav ne "verkko-base") &&
      ($flav ne "verkko-trio") &&
      ($flav ne "verkko-hi-c") &&
      ($flav ne "verkko-thic")) {
    $iflavor  = 1;
  }

  $missing = $missingCH . $missingCT . $missingB . $missingT . $missingH . $missingC;

  #  If in 'verkko-full' mode, run all verkko assemblies as a set of jobs:
  #      base -> trio -> thic
  #           -> hi-c
  if ($flav eq "verkko-full") {
    my $scrB = createVerkkoBase   ($samp, $hifi, $nano,               $missingB,        "", $compl);
    my $scrT = createVerkkoTrio   ($samp, $hifi, $nano,               $missingT, $unavailT, $compl);   #  Trio DBs implicitly exist.
    my $scrH = createVerkkoHiC    ($samp, $hifi, $nano, $hic1, $hic2, $missingH, $unavailH, $compl);
    my $scrC = createVerkkoTrioHiC($samp, $hifi, $nano, $hic1, $hic2, $missingC, $unavailC, $compl);

    my $jidB = submitIf($samp, "verkko-base", $scrB, $complB, undef, $missingB,        "", $$opts{"submit"});   #  Depends on no earlier job.
    my $jidT = submitIf($samp, "verkko-trio", $scrT, $complT, $jidB, $missingT, $unavailT, $$opts{"submit"});   #  Depends on base.
    my $jidH = submitIf($samp, "verkko-hi-c", $scrH, $complH, $jidB, $missingH, $unavailH, $$opts{"submit"});   #  Depends on base.
    my $jidC = submitIf($samp, "verkko-thic", $scrC, $complC, $jidT, $missingC, $unavailC, $$opts{"submit"});   #  Depends on base AND trio.
  }

  #  Otherwise, run a single flavor.
  else {
    my $scr;
    my $unavail;

    if ($flav eq "canu-hifi")   { $scr = createCanuHiFi     ($samp, $hifi,                      $missing, $unavail =        "", $compl); }
    if ($flav eq "canu-trio")   { $scr = createCanuTrio     ($samp,        $nano, $mati, $pati, $missing, $unavail = $unavailT, $compl); }

    if ($flav eq "verkko-base") { $scr = createVerkkoBase   ($samp, $hifi, $nano,               $missing, $unavail =        "", $compl); }
    if ($flav eq "verkko-trio") { $scr = createVerkkoTrio   ($samp, $hifi, $nano,               $missing, $unavail = $unavailT, $compl); }   #  Trio DBs implicitly exist.
    if ($flav eq "verkko-hi-c") { $scr = createVerkkoHiC    ($samp, $hifi, $nano, $hic1, $hic2, $missing, $unavail = $unavailH, $compl); }
    if ($flav eq "verkko-thic") { $scr = createVerkkoTrioHiC($samp, $hifi, $nano, $hic1, $hic2, $missing, $unavail = $unavailC, $compl); }

    if    ($compl)              { print "$samp/$flav - FINISHED\n";                                        }
    elsif (-e "$scr.jid")       { print "$samp/$flav - RUNNING\n";                                         }
    elsif (-e "$scr.err")       { print "$samp/$flav - CRASHED\n";                                         }
    elsif ($iflavor)            { print "$samp/$flav - INVALID-FLAVOR '$flav'\n";                          }
    elsif ($unavail)            { print "$samp/$flav - UNAVAILABLE ($unavail)\n";            return undef; }
    elsif ($missing)            { print "$samp/$flav - MISSING-INPUTS (missing $missing)\n";               }
    elsif (! $$opts{"submit"})  { print "$samp/$flav - READY-TO-COMPUTE\n";                                }
    else                        { print "$samp/$flav - SUBMITTED\n";  system("sbatch $scr.sh > $scr.jid"); }
  }
}

1;
