#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
#use Cwd 'abs_path';
#use Cwd;

##### ##### ##### ##### #####

use Getopt::Std;
use vars qw( $opt_a $opt_v);


# Usage
my $usage = "
fasta2nuccomp.pl - get nuccomp from fasta file.
		      by
		Brian J. Knaus
		  April 2012

Copyright (c) 2012 Brian J. Knaus.
License is hereby granted for personal, academic, and non-profit use.
Commercial users should contact the author (http://brianknaus.com).

Great effort has been taken to make this software perform its said
task however, this software comes with ABSOLUTELY NO WARRANTY,
not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

Usage: perl fasta2nuccomp.pl options
 required:
  -a	fasta file.
 optional:
  -v	verbose mode [optional T/F, default is F].

";

# command line processing.
getopts('a:v:');
die $usage unless ($opt_a);

my ($inf, $verb);

$inf	= $opt_a if $opt_a;
$verb	= $opt_v ? $opt_v : "F";

##### ##### ##### ##### #####
# Globals.

my (@names, @len, @at, @gc, @n, @poly, @other);
my $outf;

#my ($in, $outf, $out);
#my @temp;

##### ##### ##### ##### #####
# Main.

# Check file.
chk();

# Manage outfile name.
$outf = outname();
#print "$outf\n";

# Process file.
proc();

#print $#at."\n";
#print $#len."\n";

# Print to text.
ptext();

# Plot.
plot();

##### ##### ##### ##### #####
# Subroutines.

sub chk{
  my @temp;
  my $in;
  open($in, "<", $inf) or die "Can't open $inf: $!";
  chomp($temp[0] = <$in>);	# First line is id.
  chomp($temp[1] = <$in>);	# Second line is sequence.
  if($temp[0] !~ /^>/){
    print $temp[0]."\n";
    print $temp[1]."\n";
    die "error: not a fasta file\n";
  }
  close $in or die "$in: $!";
}

sub outname{
  my @temp = fileparse($inf, (".fa",".fasta"));
  return $temp[0];
#  return $temp[0]."_nuccomp.txt";
}

sub proc{
  my @temp;
  my $in;
  open($in, "<", $inf) or die "Can't open $inf: $!";
  chomp($temp[1] = <$in>);	# Initialize.
  $temp[2] = "";	# Initialize.
  while(<$in>){
    chomp($temp[0] = $_);
    if($temp[0] =~ /^>/){
    # New sequence.
#    print "$temp[1]\n";
    push @names, substr $temp[1],1;
    nuccomp($temp[2]);
    $temp[1] = $temp[0];
    $temp[2] = "";
    } else {
      # Concatenate sequence.
      $temp[2] = $temp[2].$temp[0];
    }
  }
#  print $#at."\n";
  # Last sequence.
  push @names, substr $temp[1],1;
  nuccomp($temp[2]); 
  close $in or die "$in: $!";
}

sub nuccomp{
  my $seq = uc(shift);
  my $len = length($seq);
  my $at = ($seq =~ tr/[AT]//);
  my $gc = ($seq =~ tr/[GC]//);
  my $n = ($seq =~ tr/N//);
  my $poly = ($seq =~ tr/[RYSWKMBDHV]//);
  my $other = $len - $at - $gc - $n - $poly;

#  print "$at\t";
#  print "$len\t";
#  print $at."\t".$#at."\t".$#len."\n";

  push @len, $len;
  push @at, $at;
  push @gc, $gc;
  push @n, $n;
  push @poly, $poly;
  push @other, $other;

#  print $#at."\t".$#len."\n";
}

sub ptext{
  my $out;
  my $outf = $outf."_nuccomp.txt";
  my $i;
  open($out, ">", $outf) or die "Can't open $outf: $!";
  print $out "Name,Length,AT,GC,IUPAC_polymorphic,N,Other\n";
  for($i = 0; $i<=$#len; $i++){
#    print "$names[$i]\n";
#    print $out "$names[$i],$len[$i],$at[$i],$gc[$i],$poly[$i],$n[$i],$other[$i]\n";
    print $out join(",", $names[$i], $len[$i], $at[$i], $gc[$i],
                    $poly[$i],$n[$i],$other[$i])."\n";
  }

#  print $out join(",",@len)."\n";
#  print $out join(",",@at)."\n";
#  print $out join(",",@gc)."\n";
#  print $out join(",",@poly)."\n";
#  print $out join(",",@n)."\n";
#  print $out join(",",@other)."\n";
#  close $out or die "$out: $!";
}

sub plot{
  my $inf = $outf."_nuccomp.txt";
  my $png = $outf."_nuccomp.png";
  my $outf = $outf.".r";
  my $out;
  open($out, ">", $outf) or die "Can't open $outf: $!";
  print $out "x <- read.table('".$inf."',sep=',',header=T)\n";
  print $out "lens <- sort(x[,'Length'], decreasing=T)\n";
  print $out "n50 <- lens[cumsum(as.numeric(lens)) >= sum(as.numeric(lens))/2][1]\n";
#  print $out "x <- read.table('".$inf."',sep=',',skip=1)\n";
  print $out "png('".$png."',width=800,height=800, pointsize=18)\n";
  print $out "par(mfrow=c(3,2))\n";

  print $out "hist(as.numeric(x[,2]),main='Length',xlab='')\n";
  print $out "legend('topright', c(paste( format(sum(as.numeric(x[,2])), big.mark=',')),'Nucleotides', paste('N50:', format(n50, big.mark=',')), paste('Min:', min(lens))),bty='n')\n";
#  print $out "legend('topright',c( paste(length(x[,2]),'sequences'), paste(sum(x[,2]),nucleotides)\n";

  print $out "hist(as.numeric(x[,3])/as.numeric(x[,2]),main='A/T',xlab='',xlim=c(0,1))\n";
  print $out "legend('topright',c(paste( format(sum(as.numeric(x[,3])), big.mark=',')),'Nucleotides'),bty='n')\n";

  print $out "hist(as.numeric(x[,5])/as.numeric(x[,2]),main='IUPAC polymorphic',xlab='',xlim=c(0,1))\n";
  print $out "legend('topright',c(paste( format(sum(as.numeric(x[,5])), big.mark=',')),'Nucleotides'),bty='n')\n";

  print $out "hist(as.numeric(x[,4])/as.numeric(x[,2]),main='G/C',xlab='',xlim=c(0,1))\n";
  print $out "legend('topright',c(paste( format(sum(as.numeric(x[,4])), big.mark=',')),'Nucleotides'),bty='n')\n";

  print $out "hist(as.numeric(x[,6])/as.numeric(x[,2]),main='Ns',xlab='',xlim=c(0,1))\n";
  print $out "legend('topright',c(paste( format(sum(as.numeric(x[,6])), big.mark=',')),'Nucleotides'),bty='n')\n";

  print $out "hist(as.numeric(x[,7])/as.numeric(x[,2]),main='Other',xlab='',xlim=c(0,1))\n";
  print $out "legend('topright',c(paste( format(sum(as.numeric(x[,7])), big.mark=',')),'Nucleotides'),bty='n')\n";

  print $out "dev.off()\n";
  close $out or die "$out: $!";

  my $cmd = "R CMD BATCH $outf";
  #print "$cmd\n";
  system $cmd;

  unlink $outf;
  unlink $outf.'.Rout';
}

##### ##### ##### ##### #####
# EOF.
