
###############################################################################
#
#  This file is part of a pipeline to (help) run Verkko assemblies on
#  numerous HPRC samples.
#
#  This is a 'United States Government Work', and is released in the public
#  domain.
#
##

#
#  Returns 0 is the assembly is not currently running.
#   - no .jid file or slurm reports the job-id is not currently active.
#   - assembly could be complete or crashed.
#   - if .jid file exists, it will be removed.
#
#  Returns 1 if the assemble is definitely currently running.
#   - .jid file exists and slurm positively reports the job is active.
#
#  Returns 2 if the assembly is possibly running, by we can't tell for sure.
#   - every other case; manual intervention required.
#

sub readJID ($) {
}

sub isRunning ($$) {
  my $samp = shift @_;
  my $flav = shift @_;
  my $dirn = "$rasm/$samp`";

  my $jid;      #  Submitted job ID we're looking for.
  my $jobid;    #  Job ID we found.
  my $state;    #
  my $reason;   #
  my $out;      #

  my $ret = 0;                          #  No .jid file
  my $msg = "idle";                     #  0 - NOT RUNNING

  if (-e "$dirn/$flav.jid") {
    open(JID, "< $dirn/$flav.jid") or die "Failed to open '$dirn/$flav.jid' for reading: $!\n";
    my $jid = int(<JID>);
    close(JID);

    open(JOBS, "dashboard_cli jobs --noheader --tab --fields jobid,state,state_reason,std_out --jobid $jid|") or die "Failed to open 'dashboard_cli jobs' for reading: $!\n";
    while (<JOBS>) {
      ($jobid, $state, $reason, $out) = split '\s+', $_;

      last if ($jobid eq $jid);
    }
    close(JOBS);

    if     (!defined($jobid)) {         #  No report found for expected jobid
      ($ret, $msg) = (2, "orphan");     #  2 - POSSIBLY RUNNING
    }
    elsif (($state eq "COMPLETED") ||   #  Job completed successfully or not.
           ($state eq "CANCELLED") ||   #  0 - NOT RUNNING
           ($state eq "FAILED") ||
           ($state eq "TIMEOUT")) {
      ($ret, $msg) = (0, "failed");
    }
    elsif (($state eq "CONFIGURING") || #  Job is pending.
           ($state eq "PENDING")) {     #  1 - RUNNING
      ($ret, $msg) = (1, "pending");
    }
    elsif (($state eq "RUNNING") ||     #  Confirmed running (or at least under slurm control)
           ($state eq "COMPLETING") ||  #  1 - RUNNING
           ($state eq "STOPPED") ||
           ($state eq "SUSPENDED")) {
      ($ret, $msg) = (1, "running");
    }
    else {                              #  Unrecognized
      ($ret, $msg) = (2, "unknown");    #  2 - POSSIBLY RUNNING
    }
  }

  return wantarray() ? ($ret, $jobid, $state, $msg) : $ret;
}

#
#  Return 0 if the assembly outputs do not exist OR if the flavor is invalid.
#  Return 1 if the assembly outputs do     exist; the assembly is complete.
#
sub isFinished ($$) {
  my $samp = shift @_;
  my $flav = shift @_;
  my $fini = 0;
  my $dirn = "$rasm/$samp";

  if     ($flav eq "canu-hifi") {
    $fini = 1   if (-e "$dirn/canu-hifi/asm.contigs.fasta");
  }
  elsif  ($flav eq "canu-trio") {
    $fini = 1   if ((-e "$dirn/canu-trio/asm-haplotypeMAT/asm.contigs.fasta") &&
                    (-e "$dirn/canu-trio/asm-haplotypePAT/asm.contigs.fasta"));
  }
  elsif  ($flav eq "verkko-base") {
    $fini = 1   if (-e "$dirn/$flav/contigs.fasta");
  }
  elsif (($flav eq "verkko-trio") ||
         ($flav eq "verkko-hi-c") ||
         ($flav eq "verkko-thic")) {
    $fini = 1   if (-e "$dirn/$flav/assembly.fasta") || (-e "$dirn/$flav/assembly.fasta.gz");
  } elsif (($flav eq "hifiasm-hi-c") ||
           ($flav eq "hifiasm-trio")) {
    $fini = 1   if (-e "$dirn/$flav/assembly.fasta") || (-e "$dirn/$flav/assembly.fasta.gz");
  }
  return $fini;
}

1;
