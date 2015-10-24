#!/usr/bin/env perl

use strict;
use warnings;

my $old_fasta_file = shift;
my $new_fasta_file = shift;

die "OLD_FASTA_FILE NEW_FASTA_FILE > NEW_SEQS.fsa \n" unless $old_fasta_file && -s $old_fasta_file && $new_fasta_file && -s $new_fasta_file;

my $orig_sep = $/;
$/ = ">";

open (IN,$old_fasta_file) || die $!;
my %hash_exist;
while (my $record=<IN>){
     chomp($record);
     next unless $record;
     my @lines = split("\n",$record);
     my $id = shift (@lines);
	#gi|444439379|ref|N
     if ($id=~/^gi\|(\d+)\|/){
	$hash_exist{$1}++;
     }
}
close IN;


open (IN,$new_fasta_file) || die $!;
while (my $record=<IN>){
     chomp($record);
     next unless $record;
     my @lines = split("\n",$record);
     my $id = shift (@lines);
	#gi|444439379|ref|N
     if ($id=~/^gi\|(\d+)\|/){
	print ">".$record if !$hash_exist{$1};
     }
}
close IN;
