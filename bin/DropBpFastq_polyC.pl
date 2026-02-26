#!/usr/bin/perl
use strict;
use warnings;

# Usage: perl DropBpFastq_polyC_stream.pl infile.fastq outfile.fastq

my ($infile, $outfile) = @ARGV;
die "Usage: $0 infile.fastq outfile.fastq\n" unless defined $outfile;

open(my $IN,  "<", $infile)  or die "Cannot open $infile: $!\n";
open(my $OUT, ">", $outfile) or die "Cannot open $outfile: $!\n";

while (1) {
  my $h = <$IN>;
  last unless defined $h;
  my $seq  = <$IN>;
  my $plus = <$IN>;
  my $qual = <$IN>;
  last unless defined $qual;

  chomp($h, $seq, $plus, $qual);

  # Basic sanity
  next unless $h =~ /^\@/;
  next unless $plus =~ /^\+/;
  next unless length($seq) == length($qual);

  # ---- Apply same trimming logic as your script ----

  # 1) remove polyG at end: G{3,} then 0-2 any chars, anchored at end
  my $orig_len = length($seq);
  $seq =~ s/G{3,}.{0,2}$//;
  my $new_len = length($seq);
  if ($new_len != $orig_len) {
    $qual = substr($qual, 0, $new_len);
  }

  # 2) remove polyC at start: 0-2 any chars then C{3,} at start
  $orig_len = length($seq);
  $seq =~ s/^.{0,2}C{3,}//;
  $new_len = length($seq);
  if ($new_len != $orig_len) {
    my $drop = $orig_len - $new_len;   # how many bases removed from front
    $qual = substr($qual, $drop);
  }

  # Print unless empty
  next if length($seq) == 0;

  print $OUT "$h\n$seq\n$plus\n$qual\n";
}

close $IN;
close $OUT;
