
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
@EXPORT = qw(awsToLocalPath awsToLocalInfo awsToLocalName fetchInfo getRemoteSize getLocalSize fetchData numFiles getFileMap getDownloadedFiles);

use strict;
use warnings "all";
no  warnings "uninitialized";

use hprc::samples;


#
#  Converts an s3 path to a local filesystem (relative) path, the first
#  for data files, the second for info files.
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
  $locf =~ s!^s3:----human-pangenomics--\w+--!hprc-data/$samp/!;
  $locf =~ s/.(fasta.gz|fastq.gz|sam|bam|cram)$//  if ($strip);

  return("$root/$locf");
}

sub awsToLocalInfo ($$) {
  my $samp  = shift @_;
  my $info  = shift @_;

  $info =~ s!/!--!g;
  $info =~ s!^s3:----human-pangenomics--\w+--!hprc-cache/aws/$samp/!;
  $info .= ".s3ls";

  return($info);
}

sub awsToLocalName ($$@) {
  my $samp  = shift @_;
  my $name  = shift @_;
  my $strip = shift @_;

  $name =~ s!/!--!g;    #  Replace /'s by --'s.
  $name =~ s!^.*--!!;   #  Delete everything up to the LAST --.
  $name =~ s/.(fasta.gz|fastq.gz|sam|bam|cram)$//  if ($strip);

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
    print STDERR "  Fetch AWS ls for $file.\n";
    system("mkdir -p hprc-cache/aws/$samp")    if (! -d "hprc-cache/aws/$samp");
    system("aws --no-sign-request s3 ls $file > $info 2> $info.err");
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

  system("mkdir -p hprc-data/$samp")   if (! -d "hprc-data/$samp");

  foreach my $type (sort keys %$types) {
    my $files = $samples{$samp}{$type};

    foreach my $f (@$files) {
      my $awsf = $f;                                     #  So we don't accidentally corrupt the file list.
      my $locf = awsToLocalPath($samp, $f);
      #y $awso = $f =~ s!^s3://human-pangenomics/!!/r;   #  Needed for s3api.

      if (-e $locf) {
        printf "%7s/%-9s EXISTS - %s\n", $samp, "$type:", $locf;
        next;
      }

      printf "%7s/%-9s FETCH  - %s\n", $samp, "$type:", $locf;

      #y $c = "aws --no-sign-request s3api get-object --bucket human-pangenomics --key '$awso' --range bytes=0-1048576 '$locf' > $locf.err 2>&1";
      #y $c = "seqrequester generate -min 1000 -max 10000 -sequences 100 -gaussian | gzip -1c > '$locf'";
      my $c = "aws --no-sign-request s3 cp '$awsf' '$locf' > $locf.err 2>&1";
      my $r = system($c);

      if ($r == 0) {                             #  If no error, remove the logging
        unlink "$locf.err";                      #  and continue on to the next file.
        next;
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
#  Return the number of files for a given sample and type.
#  This is NOT the number of files that are present locally;
#  use getFileMap() to get that list.
#
sub numFiles ($$) {
  my $samp  = shift @_;
  my $type  = shift @_;

  $type = "hifi"   if ($type eq "hifi-cutadapt");

  my $fles  = $samples{$samp}{$type};
  return (defined($fles)) ? scalar(@$fles) : 0
}


#
#  For a given sample/read-type, returns a map from aws-name to local-name.
#  If the file hasn't been downloaded, local-name is the empty string.
#
sub getFileMap ($$) {
  my $samp  = shift @_;
  my $type  = shift @_;
  my $hics;
  my $hici;
  my $subd;

  if      ($type eq "hic1") {       #  If asked for a specific Hi-C end, make
    $hics = "_R1_001.fastq.gz";     #  two extensions, one that we want to SAVE,
    $hici = "_R2_001.fastq.gz";     #  one that we want to INGORE.  We'll flag
    $type = "hic";                  #  anything not in either of those as 'missing',
  } elsif ($type eq "hic2") {       #
    $hics = "_R2_001.fastq.gz";     #  Also use 'hic' as the original type for
    $hici = "_R1_001.fastq.gz";     #  file discovery.
    $type = "hic";
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
      if    ($locn =~ m/$hici$/)  { next;                                  }
      elsif ($locn !~ m/$hics$/)  { die "Confused by hic name in $locn\n"; }
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

      #if (! -e $locf) {
      #  $locn =~ s/.(f(ast){0,1}[aq].gz|sam|bam|cram)$//i;   #  Make reported error name reflect
      #  $locn =~ s!/$samp/!/$samp/$subd/!;                   #  what we're searching for,
      #}
    }

    #  Now insert in the map based on if the file exists or not.

    if (-e $locf) { $localmap{$awsf} = $locf; }
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
    return ""   if (! -e $f);    #  if not, return empty string.
  }

  return join " ", values %m;
}

1;
