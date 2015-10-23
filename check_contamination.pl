#!/usr/bin/env perl 

=pod

=head1 NAME 

 check_contamination.

=head1 USAGE

 -i      =>     Genome Contig File (not scaffolds unless no other choice)
 -d      =>     Path to archae_bacteria_plasmid_viral Kraken Database

=cut

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;

my $kraken_exec  = &find_program('kraken');
my ($fasta_file,$database_dir,$number_of_seqs,%split_file_data);
my $kmer_size = 25;
my $cpus = 4;
my $split_size = 5000;

GetOptions(
	'infile=s' => \$fasta_file,
	'dir=s' => \$database_dir,
	'cpus:i' => \$cpus,
	'size:i' => \$split_size,
);

pod2usage "No genome file provided" unless $fasta_file && -s $fasta_file;
pod2usage "Invalid Kraken database option\n" unless $database_dir && -d $database_dir && -s "$database_dir/database.kdb" && -s "$database_dir/database.idx";

&split_file($fasta_file);
$number_of_seqs=`grep -hc "^>" $fasta_file.split` if !$number_of_seqs;
die "No sequences split\n" unless $number_of_seqs && $number_of_seqs > 0;
chomp($number_of_seqs);

if (!-s "check_outcomes.labels"){
	print "Processing $number_of_seqs sequences with $kraken_exec and $database_dir\n";
	my $cmd_options = " --db $database_dir --threads $cpus --fasta-input --classified-out contaminants.fasta --output check_outcomes --only-classified-output $fasta_file.split";
	system($kraken_exec.$cmd_options);
	system($kraken_exec."-filter --threshold 0.05 --db $database_dir check_outcomes > check_outcomes.0.05"); 
	system($kraken_exec."-translate  --mpa-format --db $database_dir check_outcomes.0.05 > check_outcomes.labels");
}
die "No kraken output!\n" unless -s "check_outcomes.labels";
&process_kraken("check_outcomes.labels");


##############################################
sub process_kraken(){
	my $file = shift;
	my $outfile = $file.".results";
	my %results;
	open (IN,$file);
	open (OUT,">$outfile");
	while (my $ln = <IN>){
		chomp($ln);
		my @data = split("\t",$ln);
		next unless $data[1] || $data[1] eq 'root';
		my $id = $data[0];
		if ($id=~s/\.(\d+)$//){
			my $subseq = $1;
			my $tax = $data[1];
			$results{$id}{$tax}++;
		}
	}

	foreach my $id (keys %results){
		print OUT "$id";
		foreach my $tax (keys %{$results{$id}}){
			print OUT "\t$tax"."\t".$results{$id}{$tax}."\t".sprintf("%.4f",$results{$id}{$tax}/$split_file_data{$id})."\t".$split_file_data{$id}."\n";
		}
	}
	close IN;
	close OUT;
	return $outfile;
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

sub find_program($) {
        my $prog = shift;
        $prog = `basename $prog`;
        chomp($prog);
        my $prog_path = `which $prog`
          || die("Can't find the $prog utility on your system\nPlease ensure it is installed and in your path\n");
        chomp($prog_path);
        return $prog_path;
}

