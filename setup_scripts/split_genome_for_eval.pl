#!/usr/bin/env perl

use strict;
use warnings;

open (my $fai, "<", shift) or syntax_help_and_exit();
open (my $regions, "<", shift) or syntax_help_and_exit();
my $outdir = shift or syntax_help_and_exit();
my $parts = shift or syntax_help_and_exit();
my $test_size = 1.0 / $parts;
my $train_size = 1.0 - $test_size;

my $total_size = $train_size + $test_size;

$train_size /= $total_size;
$test_size /= $total_size;


my @chrom_ids = ();
my %chrom_sizes = ();
my @core_chr_ids = ();
my @core_pos = ();
my @core_end = ();
my @core_acc_size = ();	


while (my $line = <$fai>) {
  my ($id, $len, $off, $npl) = split(/\t/, $line);
  push @chrom_ids, $id;
  $chrom_sizes{$id} = $len; 
  print STDERR "Chrom $id $len\n";
}
close $fai;

while (my $line = <$regions>) {
  my ($id, $start, $stop, $type) = split(/\t/, $line);
  next unless $type =~ /(^core)|(^cons)/i;
  push @core_chr_ids, $id;
  push @core_pos, $start;
  push @core_end, $stop;
  push @core_acc_size, $stop - $start + 1 + ($#core_acc_size >= 0 ? $core_acc_size[$#core_acc_size] : 0);
}

my $total_bps = $#core_acc_size >= 0 ? $core_acc_size[$#core_acc_size] : 0;

my $train_bps = $train_size * $total_bps;
my $test_bps = $test_size * $total_bps;

print STDERR "Total bps in core/conserverved regions is $total_bps\n";
print STDERR "Training bps target is $train_bps\n";
print STDERR "Test bps target is $test_bps\n";

my $test_prev = 0;
my $test_next = 0;
my $index = 0;
my @intervals = ();
my $start_i = 0;
my $acc_target = $test_bps;
my $last_close = 0;
my $set_count = 0;

my @set_ranges = ();
print STDERR "First pass...\n";

for my $i (0 .. $#core_chr_ids) {
  $test_next = $core_acc_size[$i];
  if ($test_next >= $acc_target) {
    if ($test_next - $acc_target <= $acc_target - $test_prev) {
      push @set_ranges, [$start_i, $i + 1, $test_next - $acc_target];
	    print STDERR "Test training set " . ++$set_count . " size is " . ($test_next - $last_close) . " " . range_str($start_i, $i + 1) . "\n";
      $last_close = $test_next;
      $start_i = $i + 1;
    } else {
      push @set_ranges, [$start_i, $i, $test_prev - $acc_target];
      print STDERR "Test training set " . ++$set_count . " size is " . ($test_prev - $last_close) . " " . range_str($start_i, $i) . "\n";
      $last_close = $test_prev;
      $start_i = $i;
    }
    
    $acc_target += $test_bps;    
  }
  $test_prev = $test_next;
}

print STDERR "Second pass...\n";
for my $i (1 .. $#set_ranges) {  
   my $left_len = $set_ranges[$i-1]->[2];
   my $right_len = $set_ranges[$i]->[2];   
   if ($left_len * $right_len < 0) {
      my $t = $set_ranges[$i-1]->[1] - 1;
      my $len = $core_end[$t] - $core_pos[$t] + 1;
      my $sum;
      if ($left_len < 0) {
        $sum = abs($left_len + $len) + abs($right_len - $len);
      } else {
        $sum = abs($left_len - $len) + abs($right_len + $len);
      }
      if ($sum < abs($left_len) + abs($right_len)) {
        if ($left_len < 0) {
	  $set_ranges[$i - 1]->[1]++;
	  $set_ranges[$i - 1]->[2] += $len;
	  $set_ranges[$i]->[0]++;
	  $set_ranges[$i]->[2] -= $len;
	} else {
	  $set_ranges[$i - 1]->[1]--;
          $set_ranges[$i - 1]->[2] -= $len;
          $set_ranges[$i]->[0]--;
          $set_ranges[$i]->[2] += $len;
	}
      }
   }
}

open (my $list, ">", "regions.list");
print $list join("\t", qw/region intervals/),"\n";

for my $i (0 .. $#set_ranges) {
    print STDERR "Test training set " . ($i + 1) . " core size is " . ($test_bps + $set_ranges[$i]->[2]) . " " . range_str($set_ranges[$i]->[0], $set_ranges[$i]->[1]) . "\n";
    open (my $out, ">" , "$outdir/test-region-$i.bed");
    my $start_i = $set_ranges[$i]->[0];
    my $end_i = $i == $#set_ranges ? $#core_chr_ids : $set_ranges[$i+1]->[0] - 1;
    my $p_start;
    my $chr;
    for my $j ($start_i .. $end_i) {
       if ($j == $start_i) {
          $chr = $core_chr_ids[$j];
	  $p_start = $j == 0 || $core_chr_ids[$j - 1] ne $chr ? 0 : int(($core_pos[$j] + $core_end[$j -1]) / 2);
       } else { 
          if ($chr ne $core_chr_ids[$j]) {
	    print $out join("\t", $chr, $p_start, $chrom_sizes{$chr}), "\n";
	    $chr = $core_chr_ids[$j];
	    $p_start = 0;
	  }
       }
       if ($j == $end_i) { 
	  my $p_end = $j == $#core_chr_ids || $core_chr_ids[$j+1] ne $chr ? $chrom_sizes{$chr} : int(($core_end[$j] + $core_pos[$j+1])/2); 
          print $out join("\t", $chr, $p_start, $p_end),"\n";
       }
    }
    close ($out);
    print join("\t", "region-$i", "test-region-$i.bed"),"\n";
}
close ($list);



sub range_str {
  my $start = shift;
  my $stop = shift;
  my $result = "";
  if ($start == 0) {
     $result = $core_chr_ids[0] . ":1 ... ";
  } elsif ($core_chr_ids[$start - 1] ne $core_chr_ids[$start]) {
     $result = $core_chr_ids[$start] . ":1 ... ";
  } else {
     $result = $core_chr_ids[$start] . ":" . int(($core_end[$start-1] + $core_pos[$start]) / 2) . " ... ";
  } 
  return $result;  
}

sub dump_intervals {
  my $index = shift;
  my $start = shift;
  my $to = shift;

  return if $to <= $start;
  
  my $spos = $start == 0 || $core_chr_ids[$start - 1] ne $core_chr_ids[$start] 

}


sub syntax_help_and_exit {
  print STDERR "split.pl fai-file regions-file train-size dest_dir number-of-parts\n";
  exit 1;
}


