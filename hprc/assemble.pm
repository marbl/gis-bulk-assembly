
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


sub submitIf ($$$$$$$) {
  my $samp    = shift @_;
  my $flav    = shift @_;
  my $scr     = shift @_;
  my $compl   = shift @_;
  my $previd  = shift @_;
  my $missing = shift @_;
  my $submit  = shift @_;
  my $jid     = undef;
  my $wait    = ($previd == 0) ? "" : "--depend=afterany:$previd";

  if    ($compl)              { print "$samp/$flav - FINISHED\n";                return undef; }
  elsif (-e "$scr.jid")       { print "$samp/$flav - RUNNING\n";                               }
  elsif (-e "$scr.err")       { print "$samp/$flav - CRASHED\n";                 return undef; }
  elsif ($missing)            { print "$samp/$flav - MISSING-INPUTS\n$missing";  return undef; }
  elsif (! $submit)           { print "$samp/$flav - READY-TO-COMPUTE\n";        return undef; }
  else                        { print "$samp/$flav - SUBMITTED\n";               system("sbatch $wait $scr.sh > $scr.jid"); }

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

  print "hifi - ", numFiles($samp, "hifi"), " $hifiMissing - $hifi\n";
  print "nano - ", numFiles($samp, "ont"),  " $nanoMissing - $nano\n";

  my $iflavor;

  if  ($flav eq "canu-hifi") {
    $missingCH  = "$samp/canu-hifi -  can't start: missing hifi-cutadapt\n"     if ($hifiMissing);
  }
  if  ($flav eq "canu-trio") {
    $missingCT  = "$samp/canu-trio -  can't start: missing ont\n"               if ($nanoMissing);
    $missingCT .= "$samp/canu-trio -  can't start: missing trio\n"              if ($trioMissing);
  }
  if (($flav eq "verkko-base") || ($flav eq "verkko-full")) {
    $missingB   = "$samp/verkko-base -  can't start: missing hifi-cutadapt\n"     if ($hifiMissing);
    $missingB  .= "$samp/verkko-base -  can't start: missing ont\n"               if ($nanoMissing);
    #missingB  .= "$samp/verkko-base -  can't start: missing trio\n"              if ($trioMissing);
    #missingB  .= "$samp/verkko-base -  can't start: missing hic\n"               if ($hicMissing);
    #missingB  .= "$samp/verkko-base -  can't start: missing base assembly\n"     if ($baseAsmMissing);
    #missingB  .= "$samp/verkko-base -  can't start: missing trio assembly\n"     if ($trioAsmMissing);
    #missingB  .= "$samp/verkko-base -  can't start: missing hapmer databases\n"  if ($hapmerMissing);
  }
  if (($flav eq "verkko-trio") || ($flav eq "verkko-full")) {
    $missingT   = "$samp/verkko-trio -  can't start: missing hifi-cutadapt\n"     if ($hifiMissing);
    $missingT  .= "$samp/verkko-trio -  can't start: missing ont\n"               if ($nanoMissing);
    $missingT  .= "$samp/verkko-trio -  can't start: missing trio\n"              if ($trioMissing);
    #missingT  .= "$samp/verkko-trio -  can't start: missing hic\n"               if ($hicMissing);
    $missingT  .= "$samp/verkko-trio -  can't start: missing base assembly\n"     if ($baseAsmMissing) && ($flav ne "verkko-full");
    #missingT  .= "$samp/verkko-trio -  can't start: missing trio assembly\n"     if ($trioAsmMissing);
    $missingT  .= "$samp/verkko-trio -  can't start: missing hapmer databases\n"  if ($hapmerMissing);
  }
  if (($flav eq "verkko-hi-c") || ($flav eq "verkko-full")) {
    $missingH   = "$samp/verkko-hi-c -  can't start: missing hifi-cutadapt\n"     if ($hifiMissing);
    $missingH  .= "$samp/verkko-hi-c -  can't start: missing ont\n"               if ($nanoMissing);
    #missingH  .= "$samp/verkko-hi-c -  can't start: missing trio\n"              if ($trioMissing);
    $missingH  .= "$samp/verkko-hi-c -  can't start: missing hic\n"               if ($hicMissing);
    $missingH  .= "$samp/verkko-hi-c -  can't start: missing base assembly\n"     if ($baseAsmMissing) && ($flav ne "verkko-full");
    #missingH  .= "$samp/verkko-hi-c -  can't start: missing trio assembly\n"     if ($trioAsmMissing);
    #missingH  .= "$samp/verkko-hi-c -  can't start: missing hapmer databases\n"  if ($hapmerMissing);
  }
  if (($flav eq "verkko-thic") || ($flav eq "verkko-full")) {
    $missingC   = "$samp/verkko-thic -  can't start: missing hifi-cutadapt\n"     if ($hifiMissing);
    $missingC  .= "$samp/verkko-thic -  can't start: missing ont\n"               if ($nanoMissing);
    $missingC  .= "$samp/verkko-thic -  can't start: missing trio\n"              if ($trioMissing);
    $missingC  .= "$samp/verkko-thic -  can't start: missing hic\n"               if ($hicMissing);
    $missingC  .= "$samp/verkko-thic -  can't start: missing base assembly\n"     if ($baseAsmMissing) && ($flav ne "verkko-full");
    $missingC  .= "$samp/verkko-thic -  can't start: missing trio assembly\n"     if ($trioAsmMissing) && ($flav ne "verkko-full");
    $missingC  .= "$samp/verkko-thic -  can't start: missing hapmer databases\n"  if ($hapmerMissing);
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
    my $scrB = createVerkkoBase   ($samp, $hifi, $nano,               $missingB, $compl);
    my $scrT = createVerkkoTrio   ($samp, $hifi, $nano,               $missingT, $compl);   #  Trio DBs implicitly exist.
    my $scrH = createVerkkoHiC    ($samp, $hifi, $nano, $hic1, $hic2, $missingH, $compl);
    my $scrC = createVerkkoTrioHiC($samp, $hifi, $nano, $hic1, $hic2, $missingC, $compl);

    my $jidB = submitIf($samp, "verkko-base", $scrB, $complB, undef, $missingB, $$opts{"submit"});   #  Depends on no earlier job.
    my $jidT = submitIf($samp, "verkko-trio", $scrT, $complT, $jidB, $missingT, $$opts{"submit"});   #  Depends on base.
    my $jidH = submitIf($samp, "verkko-hi-c", $scrH, $complH, $jidB, $missingH, $$opts{"submit"});   #  Depends on base.
    my $jidC = submitIf($samp, "verkko-thic", $scrC, $complC, $jidT, $missingC, $$opts{"submit"});   #  Depends on trio.
  }

  #  Otherwise, run a single flavor.


#
#
#  SUPPORT for non-trio hi-c only
#
#

  else {
    my $scr;

    if ($flav eq "canu-trio")   { $scr = createCanuTrio     ($samp,        $nano, $mati, $pati, $missing, $compl); }
    if ($flav eq "canu-hifi")   { $scr = createCanuHiFi     ($samp, $hifi,                      $missing, $compl); }

    if ($flav eq "verkko-base") { $scr = createVerkkoBase   ($samp, $hifi, $nano,               $missing, $compl); }
    if ($flav eq "verkko-trio") { $scr = createVerkkoTrio   ($samp, $hifi, $nano,               $missing, $compl); }   #  Trio DBs implicitly exist.
    if ($flav eq "verkko-hi-c") { $scr = createVerkkoHiC    ($samp, $hifi, $nano, $hic1, $hic2, $missing, $compl); }
    if ($flav eq "verkko-thic") { $scr = createVerkkoTrioHiC($samp, $hifi, $nano, $hic1, $hic2, $missing, $compl); }

    if    ($compl)              { print "$samp/$flav - FINISHED\n";         }
    elsif (-e "$scr.jid")       { print "$samp/$flav - RUNNING\n";          }
    elsif (-e "$scr.err")       { print "$samp/$flav - CRASHED\n";          }
    elsif ($iflavor)            { print "$samp/$flav - INVALID-FLAVOR '$flav'\n";  }
    elsif ($missing)            { print "$samp/$flav - MISSING-INPUTS\n$missing";  }
    elsif (! $$opts{"submit"})  { print "$samp/$flav - READY-TO-COMPUTE\n"; }
    else                        { print "$samp/$flav - SUBMITTED\n"; system("sbatch $scr.sh > $scr.jid"); }
  }
}

1;
