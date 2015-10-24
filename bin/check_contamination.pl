#!/usr/bin/env perl 

=pod

=head1 NAME 

 check_contamination.

=head1 USAGE

 -i      =>     Genome Contig File (not scaffolds unless no other choice)
 -d      =>     Path to archae_bacteria_plasmid_viral Kraken Database
 -cpu    =>     Number of threads for kraken (def 4)
 -split  =>     Split long sequences in this many bp + kmer size (def. 5000).
 -kmer   =>     Kmer to use. Fixed at 31.
 -prop   =>     Minimum % proportion of subsequences that needs to be positive (def. 20)

=cut

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;

my $kraken_exec  = &find_program('kraken');
my ($fasta_file,$database_dir,$number_of_seqs,%split_file_data);
my $kmer_size = 31;
my $cpus = 4;
my $split_size = 5000;
my $prop_cutoff = 20; # perc

GetOptions(
	'infile=s' => \$fasta_file,
	'dir=s' => \$database_dir,
	'cpus:i' => \$cpus,
	'size:i' => \$split_size,
	'prop_cut:i' =>\$prop_cutoff
);

pod2usage "No genome file provided" unless $fasta_file && -s $fasta_file;
pod2usage "Invalid Kraken database option\n" unless $database_dir && -d $database_dir && -s "$database_dir/database.kdb" && -s "$database_dir/database.idx";

if (-s "$fasta_file.split"){
  &check_split("$fasta_file.split");
}else{
  &split_file($fasta_file);
}

die "No sequences split\n" unless $number_of_seqs && $number_of_seqs > 0;
chomp($number_of_seqs);

if (!-s "check_outcomes.labels"){
	print "Processing $number_of_seqs sequences with $kraken_exec and $database_dir\n";
	my $cmd_options = " --db $database_dir --threads $cpus --fasta-input --output check_outcomes --only-classified-output $fasta_file.split";
	system($kraken_exec.$cmd_options);
	system($kraken_exec."-filter --threshold 0.05 --db $database_dir check_outcomes > check_outcomes.0.05"); 
	system($kraken_exec."-translate  --mpa-format --db $database_dir check_outcomes.0.05 > check_outcomes.labels");
}
die "No kraken output!\n" unless -s "check_outcomes.labels";
&process_kraken("check_outcomes.labels");


##############################################
sub check_split(){
	%split_file_data=();
	my $file = shift;
	open (IN,$file) || die $!;
	my $orig_sep = $/;
	$/ = ">";
	while (my $record=<IN>){
		chomp($record);
		next unless $record;
		my @lines = split("\n",$record);
		my $id = shift (@lines);
		$id=~s/\s+.+$//;
		if ($id=~s/\.(\d+)$//){
			my $counter = $1;
			$number_of_seqs++;
			$split_file_data{$id}++;
		}else{die;}
	}
	close IN;
	$/ = $orig_sep;
	my $seq_count = scalar(keys %split_file_data);
	print "Found: $number_of_seqs from $seq_count            \n";
}


sub split_file(){
	my ($file) = @_;
	open (IN,$file) || die $!;
	open (OUT,">$file.split");
	print "Splitting files\n";
	my $orig_sep = $/;
	$/ = ">";
	$|=1;
	my $seq_count = 0;
	while (my $record=<IN>){
		chomp($record);
		next unless $record;
		$seq_count++;
		my @lines = split("\n",$record);
		my $id = shift (@lines);
		$id =~s/\s+.+$//;
		my $seq = join('',@lines);
		$seq=~s/\s+//g;
		$seq=~s/[^ATGCN]/N/g;
		my $seq_length = length($seq);
		my @seq_array = split('',$seq);
		my $counter = 1;
		for (my $i=0;$i<$seq_length;$i+=$split_size){
			last if $i >= $seq_length;
			my $k = $i + $split_size - 1 + $kmer_size;
			$k = $seq_length-1 if $k > $seq_length-1;
			my $subseq = join('',@seq_array[$i..$k]);
			print OUT ">$id.$counter ".($i+1).':'.($k+1)."\n$subseq\n";
			$counter++;
			$number_of_seqs++;
			$split_file_data{$id}++;
			last if $k >= $seq_length;
		}
		print "Created: $number_of_seqs from $seq_count contigs           \r" if $seq_count % 10 == 0;
	} 
	close IN;
	close OUT;
	$/=$orig_sep;
	$|=0;
	print "Created: $number_of_seqs from $seq_count            \n";
}

sub process_kraken(){
	my $file = shift;
	my $outfile = $file.".results";
	my %results;
	open (IN,$file);
	while (my $ln = <IN>){
		chomp($ln);
		my @data = split("\t",$ln);
		next unless $data[1] || $data[1] eq 'root';
		my $id = $data[0];
		if ($id=~s/\.(\d+)$//){
			my $subseq = $1;
			my $tax = $data[1];
			next if $tax eq 'root';
			$results{$id}{$tax}++;
		}
	}
	open (OUT,">$outfile");
	print OUT "ID\tTAXONOMY\tSUBSEQ_TAX_HITS\tSUBSEQ_PROPORTION_%\tTOTAL_SUBSEQS\n";
	foreach my $id (keys %results){
		foreach my $tax (keys %{$results{$id}}){
			my $prop = sprintf("%.4f",$results{$id}{$tax}/$split_file_data{$id})*100;
			print OUT $id."\t$tax"."\t".$results{$id}{$tax}."\t$prop\t".$split_file_data{$id}."\n" if $prop > $prop_cutoff;
		}
	}
	close IN;
	close OUT;
	return $outfile;
}



sub find_program($) {
        my $prog = shift;
        $prog = `basename $prog`;
        chomp($prog);
        my $prog_path = `which $prog`
          || die("Can't find the $prog utility on your system\nPlease ensure it is installed and in your path\n");
        chomp($prog_path);
        return $prog_path;
}

