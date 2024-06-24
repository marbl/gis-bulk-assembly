
###############################################################################
#
#  This file is part of a pipeline to (help) run Verkko assemblies on
#  numerous HPRC samples.
#
#  This is a 'United States Government Work', and is released in the public
#  domain.
#
##

package hprc::status;
require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(checkStatus reportSlurmTime);

use strict;
use warnings "all";
no  warnings "uninitialized";

use Time::Local;
use Date::Parse;

use hprc::samples;
use hprc::assemble;


sub checkStatus ($$) {
  my $samp = shift @_;
  my $opts = shift @_;
  my $flav = $$opts{"flavor"};
  my $dirn = "$rasm/$samp";
  my @logs;
  my $showJobs  = exists($$opts{'jobs'});    #  What stuff to show.
  my $showTimes = exists($$opts{'times'});

  #  Report simple status unless 'jobs' or 'times' is explicitly requested.

  if (!$showJobs && !$showTimes) {
    if    (isFinished($samp, $flav)) {
      print "$samp/$flav\tFINISHED\n";
    }
    elsif (-e "$dirn/$flav.jid") {
      print "$samp/$flav\tRUNNING\n";
    }
    elsif (-e "$dirn/$flav.sh") {
      print "$samp/$flav\tFAILED (or not submitted)\n";
    }
    else {
      print "$samp/$flav\tNOT STARTED\n";
    }

    return;
  }      

  #  Find logs.

  open(FIL, "find $rasm -maxdepth 2 -type f |");
  while (<FIL>) {
    s/^\s+//;
    s/\s+$//;

    if (m!$samp/$flav\.\d+\.log$!) {
      push @logs, "$_";
    }
  }
  close(FIL);

  #  Nothing?

  return      if (scalar(@logs) == 0);

  #  Get slurm info without calling sacct over and over.
  #    - scan logs to find the earliest/latest job dates.
  #    - request info on all jobs for those dates.
  #
  my $bgnTime = 1999999999;   my $bgnDate = "Tue May 17 23:33:20 2033";
  my $endTime = 1000000000;   my $endDate = "Sat Sep  8 21:46:40 2001";

  foreach my $log (@logs) {
    #print STDERR "Scan $log for dates.\n";

    open(LOG, "< $log") or die "Failed to open log '$log' for reading: $!\n";
    while (<LOG>) {
      chomp;

      if (m/202\d\]$/) {  #  This is so we can catch errors parsing dates.
        if (m/^\s*\[(Sun|Mon|Tue|Wed|Th[ur]|Fri|Sat)\s(...\s+\d+\s\d\d:\d\d:\d\d\s\d\d\d\d)\]$/) {
          my  $tt = str2time($2);
          my ($ss, $mm, $hh, $da, $mo, $yr, $zn, $cn) = strptime($2);

          if ($tt < $bgnTime)  { $bgnTime = $tt;  $bgnDate = sprintf "%4d-%02d-%02dT%02d:%02d:%02d", $yr+1900, $mo+1, $da,  0,  0,  0;  }
          if ($endTime < $tt)  { $endTime = $tt;  $endDate = sprintf "%4d-%02d-%02dT%02d:%02d:%02d", $yr+1900, $mo+1, $da, 23, 59, 59;  }
        }
        else {
          die "Failed to parse date: '$_'\n";
        }
      }
    }
    close(LOG);
  }

  #print "$bgnTime - $bgnDate\n";
  #print "$endTime - $endDate\n";

  my %jobToStart;
  my %jobToEnd;
  my %jobToQTime;
  my %jobToElapsed;
  my %jobToMaxMem;
  my %jobToAvgCPU;

  open(SLM, "dashboard_cli jobs --noheader --null 0 --archive --since $bgnDate --until $endDate --raw --fields jobid,start_time,end_time,queued_time,elapsed_time,mem_max,cpu_avg |");
  while (<SLM>) {
    chomp; 
    if (m/^(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+\.?\d*)$/) {
      printf "%10s -> start=%10s end=%10s queue=%10s elapsed=%10s maxmem=%10s avgcpu=%10s\n", $1, $2, $3, $4, $5, $6, $7;
      $jobToStart{$1}   = $2;
      $jobToEnd{$1}     = $3;
      $jobToQTime{$1}   = $4;
      $jobToElapsed{$1} = $5;
      $jobToMaxMem{$1}  = $6;
      $jobToAvgCPU{$1}  = $7;
    } else {
      print STDERR "WARNING: failed to parse dashboard_cli line '$_'\n";
    }
  }
  close(SLM);

  #  Scan logs, for real this time.

  if ($showJobs == 1) {
    printf "%-42s" .                                 "  ------------ SUBMIT ------------  ----------- RESUBMIT -----------  ---------- COMPLETION ----------\n", "$samp/$flav";
    printf "snakeID     rule name                       slurmID                     date  slurmID                     date  slurmID                     date\n";
    printf "----------  ------------------------------  --------------------------------  --------------------------------  --------------------------------\n";
  }

  my %stageCPU;
  my %stageWmin;
  my %stageWmax;
  my %stageWall;
  my %stageQueue;
  my %stageMaxMem;

  foreach my $log (@logs) {
    my $version;
    my $date;
    my $rule;
    my $jobnum;

    my $submit;

    my @tids;

    my %tid2jid;
    my %tid2rule;

    my %submitted;
    my %resubmit;
    my %completed;

    open(LOG, "< $log") or die "Failed to open log '$log' for reading: $!\n";
    while (<LOG>) {
      chomp;

      if (m/^Launching\sverkko\s(.*)$/) {
        $version = $1;
      }

      if (m/^\[(Sun|Mon|Tue|Wed|Thr|Fri|Sat)\s(...)\s+(\d+)\s(\d\d:\d\d:\d\d)\s(\d\d\d\d)\]$/) {
        $date = sprintf("%4d-%03s-%02d %s", $5, $2 =~ tr/a-z/A-Z/r, $3, $4);
      }

      if ((m/^rule\s+(.*):/) ||
          (m/^localrule\s+(.*):/) ||
          (m/^checkpoint\s+(.*):/)) {
        $rule = $1;
      }

      #
      #  Job completion.
      #

      if ((m/^Finished\sjob\s(\d+).$/) ||
          (m/^Error\sexecuting.*jobid:\s+(\d+),\s/)) {
        my $tid = $1;  #sprintf("tid%04d", $1);
        my $jid = $tid2jid{$tid};
        my $jid1;
        my $jid2;
        my $status;

        if (m/Finished/)  {  $status = "OK  "; }
        if (m/Error/)     {  $status = "FAIL"; }

        if ($jid =~ m/^([jid]{0,3}\d+),([jid]{0,3}\d+)$/) {
          $jid1 = $1;
          $jid2 = $2;

          $completed{$tid} = $date;
        }
        else {
          $jid1 = $jid;
          $jid2 = $jid;

          $completed{$tid} = $date;
        }

        my $s = $submitted{$tid};
        my $r = $resubmit {$tid};
        my $c = $completed{$tid};

        if ($jid1 ne "(local)") {
          my $r = $tid2rule{$tid};

          $stageCPU   {$r} += ($jobToElapsed{$jid1} * $jobToAvgCPU{$jid1});
          $stageWmin  {$r}  = $jobToStart{$jid1}     if (!exists($stageWmin{$r})) || ($stageWmin{$r} > $jobToStart{$jid1});
          $stageWmax  {$r}  = $jobToEnd{$jid1}       if (!exists($stageWmax{$r})) || ($stageWmax{$r} < $jobToEnd{$jid1});
          $stageWall  {$r} += $jobToElapsed{$jid1};
          $stageQueue {$r} += $jobToQTime{$jid1};
          $stageMaxMem{$r}  = $jobToMaxMem{$jid1}    if ($stageMaxMem{$r} < $jobToMaxMem{$jid1});
        }

        if (($jid1 ne $jid2) && ($jid1 ne "(local)")) {
          my $r = $tid2rule{$tid};

          $stageCPU   {$r} += ($jobToElapsed{$jid2} * $jobToAvgCPU{$jid2});
          $stageWmin  {$r}  = $jobToStart{$jid2}     if (!exists($stageWmin{$r})) || ($stageWmin{$r} > $jobToStart{$jid2});
          $stageWmax  {$r}  = $jobToEnd{$jid2}       if (!exists($stageWmax{$r})) || ($stageWmax{$r} < $jobToEnd{$jid2});
          $stageWall  {$r} += $jobToElapsed{$jid2};
          $stageQueue {$r} += $jobToQTime{$jid2};
          $stageMaxMem{$r}  = $jobToMaxMem{$jid2}    if ($stageMaxMem{$r} < $jobToMaxMem{$jid2});
        }

        if ($showJobs == 1) {
          if (defined($r)) {
            printf "%-5s %s  %-30s  %-11s %-20s  %-11s %-20s  %-11s %-20s\n",
                $tid,
                $status,
                $tid2rule{$tid},
                $jid1, $s,
                $jid2, $r,
                $jid2, $c;
          }
          else {
            printf "%-5s %s  %-30s  %-11s %-20s  %-11s %-20s  %-11s %-20s\n",
                $tid,
                $status,
                $tid2rule{$tid},
                $jid1, $s,
                "", "",
                $jid2, $c;
          }
        }

        delete $tid2jid{$tid};
        delete $tid2rule{$tid};
        delete $submitted{$tid};
        delete $resubmit {$tid};
        delete $completed{$tid};
      }

      #
      #  Job submission or resubmission.
      #

      if ((m/^localrule/) ||
          (m/^localcheckpoint/)) {
        while ($_ !~ m/^\s+jobid:\s+(\d+)$/) {
            $_ = <LOG>;
        }

        if ($_ =~ m/^\s+jobid:\s+(\d+)$/) {
          my $tid = $1;  #sprintf("tid%04d", $1);

          $tid2jid{$tid}   = "(local)";
          $tid2rule{$tid}  = $rule;
          $submitted{$tid} = $date;
        }
      }

      if (m/^Submitted\sjob\s(\d+)\swith\sexternal\sjobid\s'(\d+)'.$/) {
        my $tid = $1;  #sprintf("tid%04d", $1);
        my $jid = $2;  #sprintf("jid%08d", $2);

        if (exists($tid2jid{$tid})) {
          $tid2jid{$tid}   = "$tid2jid{$tid},$jid";
          $resubmit{$tid}  = $date;
        }
        else {
          $tid2jid{$tid}   = $jid;
          $tid2rule{$tid}  = $rule;
          $submitted{$tid} = $date;
        }
      }
    }
    close(LOG);


    foreach my $tid (keys %tid2jid) {
      my $jid = $tid2jid{$tid};
      my $jid1;
      my $jid2;

      if ($jid =~ m/^([jid]{0,3}\d+),([jid]{0,3}\d+)$/) {
        $jid1 = $1;
        $jid2 = $2;
      }
      else {
        $jid1 = $jid;
        $jid2 = $jid;
      }

      my $s = $submitted{$tid};
      my $r = $resubmit {$tid};
      my $c = $completed{$tid};

      if ($showJobs == 1) {
        if (defined($r)) {
          printf "%-5s RUN   %-30s  %-11s %-20s  %-11s %-20s  %-11s %-20s\n",
              $tid,
              $tid2rule{$tid},
              $jid1, $s,
              $jid2, $r,
              "", "";
        }
        else {
          printf "%-5s RUN   %-30s  %-11s %-20s  %-11s %-20s  %-11s %-20s\n",
              $tid,
              $tid2rule{$tid},
              $jid1, $s,
              "", "",
              "", "";
        }
      }
    }

    print "\n";
    print "\n";
  }

  if ($showTimes == 1) {
    printf "%-30s" .                     "       Actual   Sequential   Sequential   Est. Total  Max Memory\n", "$samp/$flav";
    printf "rule name                          Wall (h)     Wall (h)    Queue (h)      CPU (h)        (GB)\n";
    printf "------------------------------ --incorrect- ------------ ------------ ------------ -----------\n";

    foreach my $r (keys %stageCPU) {
      printf "%-30s %12.2f %12.2f %12.2f %12.2f %11.3f\n", $r, 
          $stageWmax  {$r} / 3600.0 - $stageWmin{$r} / 3600.0,
          $stageWall  {$r} / 3600.0,
          $stageQueue {$r} / 3600.0,
          $stageCPU   {$r} / 3600.0,
          $stageMaxMem{$r} / 1024.0;
    }
  }
}

1;
