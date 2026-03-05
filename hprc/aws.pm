
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
@EXPORT = qw(awsToLocalPath awsToLocalInfo awsToLocalName fetchInfo getRemoteSize getLocalSize fetchData getFileMap getDownloadedFiles getCoverage);

use strict;
use warnings "all";
no  warnings "uninitialized";

use hprc::samples;


#
#  Converts an s3 path to a local filesystem (relative) path, the first
#  for data files, the second for info files. Also handles ftp paths from
#  the sra and paths on the local server but outside the assembly workspace.
#    Data files are in 'hprc-data/$samp/'.
#    Info files are in 'hprc-cache/aws/$samp/' and have '.s3ls' appended.
#
#  awstoLocalName() strips off all the aws and encoded directory paths
#  and returns just the name of the file.
#
sub awsToLocalPath ($$@) {
  my $samp  = shift @_;
  my $locf  = shift @_;
  my $strip = shift @_;

  $locf =~ s!/!--!g;

  if ($locf =~ /^s3:/) {
    $locf =~ s!^s3:----human-pangenomics--\w+--!!;
    # for s3 paths not from human-pangenomics bucket
    $locf =~ s!^s3:----!!;
  }
  elsif ($locf =~ /^ftp:/) {
    # for ftp links from sra.ebi.ac.uk
    $locf =~ s!^ftp:----ftp.sra.ebi.ac.uk--!!;
    #for ftp links from 1000genomes
    $locf =~ s!^ftp:----ftp.1000genomes.ebi.ac.uk--!!;
  }
  elsif ($locf =~ /^local:/) {
    $locf =~ s!^.*--!!;
  }

  $locf = (length($locf)) > 255-10 ? substr($locf, length($locf)-255+10, 255) : $locf;
  $locf = "$data/$samp/$locf";
  $locf =~ s/.(fasta.gz|fastq.gz|sam|bam|cram|fq.gz|fa.gz)$//  if ($strip);
  die "Error: path $locf is too long" if length($locf) > 4096;

  return("$locf");
}

sub awsToLocalInfo ($$) {
  my $samp  = shift @_;
  my $info  = shift @_;

  $info =~ s!/!--!g;
  if ($info =~ /^s3:/) {
    $info =~ s!^s3:----human-pangenomics--\w+--!!;
    # for s3 paths not from human-pangenomics bucket
    $info =~ s!^s3:----!!;
  }
  elsif ($info =~ /^ftp:/) {
    # for ftp links from sra.ebi.ac.uk
    $info =~ s!^ftp:----ftp.sra.ebi.ac.uk--!!;
    # for ftp links from 1000genomes
    $info =~ s!^ftp:----ftp.1000genomes.ebi.ac.uk--!!;
  }
  elsif ($info =~ /^local:/) {
    $info =~ s!^.*--!!;
  }

  # leave space to append .s3ls.err to the filename
  $info = (length($info))     > 255-10 ? substr($info, 0, 255-10) : $info;
  $info = "hprc-cache/aws/$samp/$info";
  $info .= ".s3ls";
  die "Error: path $info is too long" if length($info) > 4096;

  return($info);
}

sub awsToLocalName ($$@) {
  my $samp  = shift @_;
  my $name  = shift @_;
  my $strip = shift @_;

  $name =~ s!/!--!g;    #  Replace /'s by --'s.
  $name =~ s!^.*--!!;   #  Delete everything up to the LAST --.
  die "Error: filename path $name is too long" if length($name) > 255;
  $name =~ s/.(fasta.gz|fastq.gz|sam|bam|cram|fq.gz|fa.gz)$//  if ($strip);

  return($name);
}


#
#  Given a sample name and a file:
#   - Fetch size info on it from aws, return the name of the info file.
#       Caches directory listings locally in 'hprc-cache/aws'.
#   - Return the aws (remote) or downloaded (local) size of the file.
#       Files not downloaded return size of zero.
#
sub fetchInfo ($$) {
  my $samp = shift @_;
  my $file = shift @_;
  my $info = awsToLocalInfo($samp, $file);

  if (! -e "$info") {
    if (($file =~ /^s3/) && system("aws --version > /dev/null 2>&1") != 0) {
      die "aws failed to run; probably need 'module load aws'.\n";
    }
    system("mkdir -p hprc-cache/aws/$samp")    if (! -d "hprc-cache/aws/$samp");
    if ($file =~ /^s3/) {
      system("aws --no-sign-request s3 ls $file > $info 2> $info.err");
    } elsif ($file =~ /^ftp/) {
      my $curl_output=`curl -I $file 2> $info.err`;
      # Parse the curl output for Last-Modified and Content-Length
      open my $o, ">", $info;
      my ($last_modified) = $curl_output =~ /Last-Modified:\s*(.*?)\s*\r?\n/;
      my ($size) = $curl_output =~ /Content-Length:\s*(\d+)\s*\r?\n/;

      # report failure if we can't get info, otherwise format it nicely
      #die "Couldn't get file info for $file, check $info.err " if (! $last_modified);
      my $formatted_date = `date -d "$last_modified" "+%Y-%m-%d %H:%M:%S"`;
      chomp($formatted_date);
      print $o "$formatted_date\t$size\t$file\n";
      close $o;
    } elsif ($file =~ /^local/) {
      my $localfile = $file;
      $localfile =~ s:^local\:/+:/:;
      # use ls to get last modified date and size
      my $ls_output=`ls -l --full-time $file 2> $info.err`;
      my @ls_fields = split /\s/, $ls_output;
      my $size = $ls_fields[4];
      my $fulldate = "$ls_fields[5] $ls_fields[6]";
      $fulldate =~ s/\..*//;
      open my $o, ">", $info;
      print $o "$fulldate\t$size\t$file\n";
      close $o;
    } else {
       die "Unknown file type $file\n"
    }
  }

  die "didn't find '$info'\n"   if (! -e "$info");

  return $info;
}

sub getRemoteSize ($$) {
  my $samp = shift @_;
  my $file = shift @_;
  my $info = fetchInfo($samp, $file);
  my $size = 0;

  open(INF, "< $info") or die "Failed to open '$info' for reading: $!\n";
  while (<INF>) {
    s/^\s+//;
    s/\s+$//;

    my @v = split '\s+', $_;

    $size += $v[2];
  }
  close(INF);

  return $size / 1024.0 / 1024.0 / 1024.0;
}

sub getLocalSize ($$) {
  my $samp = shift @_;
  my $file = shift @_;
  my $size;

  if (-e $file)   { $size = -s $file                         }   #  $file is local file
  else            { $size = -s awsToLocalPath($samp, $file); }   #  $file is aws-file path

  return $size / 1024.0 / 1024.0 / 1024.0;
}



#
#  Given a sample name and a file type, fetches all the files into
#  the local filesystem.
#
sub fetchData ($$$) {
  my $samp  = shift @_;
  my $types = shift @_;
  my $opts  = shift @_;

  system("mkdir -p $data/$samp")   if (! -d "$data/$samp");

  foreach my $type (sort keys %$types) {
    my %fileMap = getFileMap($samp, $type, "even those that don't exist");
    my @files = keys %fileMap; #$samples{$samp}{$type};

    next  if (!(@files));

    foreach my $f (sort @files) {
      # "AWS" name could also be an ftp path ("ftp://") or a local path ("local://")
      my $awsf = $f;                                     #  So we don't accidentally corrupt the file list.
      my $locf = awsToLocalPath($samp, $f);
      #y $awso = $f =~ s!^s3://human-pangenomics/!!/r;   #  Needed for s3api.

      if (-e $locf) {
        printf "%7s/%-9s EXISTS - %s\n", $samp, "$type:", $locf;
        next;
      }
     # also check for cutadapt version
     my $loccut = awsToLocalPath("$samp/hifi-cutadapt", $f, 1) . ".fasta.gz";
     if (-e $loccut ) {
        printf "%7s/%-9s EXISTS CUTADAPT - %s\n", $samp, "$type:", $loccut;
        next;
     }

      printf "%7s/%-9s FETCH  - %s\n", $samp, "$type:", $locf;

      if (($awsf =~ /^s3/) && (system("aws --version > /dev/null 2>&1") != 0)) {
        die "aws failed to run; probably need 'module load aws'.\n";
      }

      #y $c = "aws --no-sign-request s3api get-object --bucket human-pangenomics --key '$awso' --range bytes=0-1048576 '$locf' > $locf.err 2>&1";
      #y $c = "seqrequester generate -min 1000 -max 10000 -sequences 100 -gaussian | gzip -1c > '$locf'";
      #my $c = "aws --no-sign-request s3 cp '$awsf' '$locf' > $locf.err 2>&1";
      my $c;
      if ($awsf =~ /^s3/) {
        $c = "aws --no-sign-request s3 cp --only-show-errors '$awsf' '$locf' > $locf.err 2>&1";
      } elsif ($awsf =~ /^ftp/) {
        $c = "curl -f -o '$locf' '$awsf' > $locf.err 2>&1";
      } elsif ($awsf =~ /^local/) {
        $awsf =~ s:^local\:/+:/:;
        $c = "ln -sf '$awsf' '$locf' > $locf.err 2>&1";
      } else {
        print "Download requested unknown file type $awsf\n";
        exit(1);
      }
      my $r = system($c);
      #print "$c\n";

      if ($r == 0) {                             #  If no error, remove the logging
        unlink "$locf.err";                      #  and continue on to the next file.
        next;
      } else {
        print "Download failed with return code $r\n";
      }

      rename '$locf',     '$locf.FAILED';        #  If an error, save the failed file and logging,
      rename '$locf.err', '$locf.FAILED.err';    #  then spit out an error and die.

      print "                  FAILED: '$c'\n";
      open(FAIL, "< $locf.err");
      $_ = <FAIL>;
      $_ = <FAIL>;
      while (<FAIL>) {
        print "                     $_";
        last;
      }
      close(FAIL);

      print "\n";

      exit(1);
    }
  }
}

#
#  For a given sample/read-type, returns a map from aws-name to local-name.
#  If the file hasn't been downloaded, local-name is the empty string.
#
sub getFileMap ($$@) {
  my $samp  = shift @_;
  my $type  = shift @_;
  my $nonE  = shift @_;   #  return even if the file doesn't exist.
  my $hics;
  my $hici;
  my $ontr10;
  my $subd;

  if      ($type eq "hic1") {       #  If asked for a specific Hi-C end, make
    $hics = qr/_L\d_1.fq.gz$|_R1_001.*fastq.gz$|001.R1.fastq.gz$|_1.fastq.gz$/;     #  two extensions, one that we want to SAVE,
    $hici = qr/_L\d_2.fq.gz$|_R2_001.*fastq.gz$|001.R2.fastq.gz$|_2.fastq.gz$/;     #  one that we want to INGORE.  We'll flag
    $type = "hic";                  #  anything not in either of those as 'missing',
  } elsif ($type eq "hic2") {       #
    $hics = qr/_L\d_2.fq.gz$|_R2_001.*fastq.gz$|001.R2.fastq.gz$|_2.fastq.gz$/;     #  Also use 'hic' as the original type for
    $hici = qr/_L\d_1.fq.gz$|_R1_001.*fastq.gz$|001.R1.fastq.gz$|_1.fastq.gz$/;     #  file discovery.
    $type = "hic";
  }
  if      ($type eq "ont-r10") {  # select only those read which are R10
     $ontr10 = qr/_R1041_|_r10.4.1_/;
     $type = "ont"
  }

  if ($type eq "hifi-cutadapt") {   #  For processed data, look in a specific
    $subd = "hifi-cutadapt";        #  subdirectory for the data file, but use
    $type = "hifi";                 #  the original type to discover files.
  }

  my $awspaths = $samples{$samp}{$type};   #  <- 'original type' used here.
  my %localmap;

  foreach my $awsf (@$awspaths) {
    my $locf = awsToLocalPath($samp, $awsf);   #  Path to local file.
    my $locn = awsToLocalName($samp, $awsf);   #  Name of local file.

    #  If returning a specific hic end, decide if the file is one we
    #  care about and present/missing or if it is one we don't care about.

    if (defined($hics)) {
      if    ($locn =~ $hici)  { next;                                  }
      elsif ($locn !~ $hics)  { die "Confused by hic name in $locn\n"; }
    }
    if (defined($ontr10)) {
      if    ($locn !~ $ontr10) { next; }
    }

    #  If a subdirectory is specified, insert it in the path AND use whatever
    #  extension exists in that directory.

    if (defined($subd)) {
      #y $of = $locf;
      my $nf = $locf;

      $nf   =~ s/.(f(ast){0,1}[aq].gz|sam|bam|cram)$//i;  #  Replace existing extension
      $nf   =~ s!/$samp/!/$samp/$subd/!;                  #  Insert new subdirectory
      $locf = undef;

      foreach my $ext (qw(fasta.gz fa.gz fastq.gz fq.gz cram bam sam)) {
        if (-e "$nf.$ext") {                              #  Find a file with the
          $locf = "$nf.$ext";                             #  same name but different
          last;                                           #  extension.
        }
      }
    }

    #  Now insert in the map based on if the file exists or not.
    if (-e $locf) { $localmap{$awsf} = $locf; }
    elsif ($nonE) { $localmap{$awsf} = $locf; }
    else          { $localmap{$awsf} = "";    }
  }

  #  And return the whole map.
  return %localmap;
}


sub getDownloadedFiles ($$) {
  my $samp  = shift @_;
  my $type  = shift @_;
  my @files;
  my $incomplete = 0;

  my %m = getFileMap($samp, $type);

  foreach my $f (values %m) {    #  Check that each files exists;
    return ""   if (! -e $f || -e "$f.err");    #  if not, return empty string.
  }

  return join " ", sort values %m;
}

sub getCoverage($$) {
  my $samp   = shift @_;
  my $types  = shift @_;

  #  Read total bases
  my $nBases = 0;

  foreach my $t (sort keys %$types) {
    my %files = getFileMap($samp, $t, "even-those-that-don't-exist");
    foreach my $file (values %files) {
      my $outs = $file;    $outs =~ s/.(fasta.gz|fastq.gz|sam|bam|cram|fq.gz|fa.gz)$/.summary/;

      if (-e "$outs") {
        open(SUM, "< $outs") or die "Failed to open '$outs': $!\n";
        while (<SUM>) {
          if (m/^00100\s+(\d+)\s+(\d+)\s+(\d+)\s+/) {
            $nBases += $3;
          }
        }
        close(SUM);
      } else {
        die "Failed to get coverage, missing summary for '$outs'. Please run read-stats first!\n";
      }
    }
  }
  return int($nBases / 3100000000);
}

1;
