
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
@EXPORT = qw(checkStatus);

use strict;
use warnings "all";
no  warnings "uninitialized";

use hprc::samples;

my $root = -e "/data/walenzbp/hprc" ? "/data/walenzbp/hprc" : ".";


#  Returns the Slurm JobID of the last known running verkko command,
#  or undef if no known last job.
#
sub findJobID ($) {
  my $samp = shift @_;
  my $jid  = undef;

  print "Query $root/assemblies/$samp.jid.\n";
  if (-e "$root/assemblies/$samp.jid") {
    open(JID, "< $root/assemblies/$samp.jid") or die "Failed to open '$root/assemblies/$samp.jid' for reading: $!\n";
    $jid = int(<JID>);
    close(JID);
  }

  die "No job id.\n" if (!defined($jid));
  return($jid);
}


#  Returns 'status' of the last known verkko command.
#  This came from verkko/src/profiles/slurm-sge-status.sh.
#
sub findJobStatus ($$) {
  my $samp = shift @_;
  my $jid  = shift @_;

  my $state = "UNKNWON";

  if (-x "/usr/local/bin/dashboard_cli") {
    my ($jobid, $state, $reason, $out);

    open(JOBS, "dashboard_cli jobs --noheader --tab --fields jobid,state,state_reason,std_out --jobid $jid|") or die "Failed to open 'dashboard_cli jobs' for reading: $!\n";
    while (<JOBS>) {
      ($jobid, $state, $reason, $out) = split '\s+', $_;

      last if ($jobid eq $jid);
    }
    close(JOBS);

    if (!defined($jobid)) {
      print STDERR "Didn't get a report for job $jid from Slurm.  Can't check if I can safely submit.\n";
      exit(0);
    }

    if    (($state eq "COMPLETED")) {
      $state = "COMPLETED";
    }
    elsif (($state eq "PENDING") ||
           ($state eq "CONFIGURING") ||
           ($state eq "RUNNING") ||
           ($state eq "SUSPENDED") ||
           ($state eq "PREEMPTED") ||
           ($state eq "COMPLETING")) {
      $state = "RUNNING";
    }
    elsif (($state eq "BOOT_FAIL") ||
           ($state eq "CANCELLED") ||
           ($state eq "DEADLINE") ||
           ($state eq "FAILED") ||
           ($state eq "NODE_FAIL") ||
           ($state eq "OUT_OF_MEMORY") ||
           ($state eq "PREEMPTED") ||
           ($state eq "TIMEOUT")) {
      $state = "FAILED";
    }
    else {
      $state = "UNKNOWN";
    }
  }

  return($state);
}


sub checkStatus ($$) {
  my $samp = shift @_;
  my $opts = shift @_;
  my @logs;

  #  Find jobid.

  #my $jid  = findJobID($samp);
  #my $stat = findJobStatus($samp, $jid);

  #  Find logs.

  open(FIL, "ls $root/assemblies/ |");   #  Note: 'ls assemblies/*log' fails if no logs present!
  while (<FIL>) {
    s/^\s+//;
    s/\s+$//;

    if (m/^$samp\.\d+\.log$/) {
      push @logs, "$root/assemblies/$_";
    }
  }
  close(FIL);

  #  Nothing?

  return      if (scalar(@logs) == 0);

  #  Scan logs.

  print "                                           ------------ SUBMIT ------------  ----------- RESUBMIT -----------  ---------- COMPLETION ----------\n";
  print "snakeID    rule name                       slurmID                     date  slurmID                     date  slurmID                     date\n";
  print "---------  ------------------------------  --------------------------------  --------------------------------  --------------------------------\n";

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

    print "\n";
    print "Scan $log\n";
    print "\n";

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
}

1;
