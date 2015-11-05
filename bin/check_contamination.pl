#!/usr/bin/env perl 

=pod

=head1 NAME 

 check_contamination.

=head1 USAGE

 -i     :s  =>     Genome Contig File (not scaffolds unless no other choice)
 -d     :s  =>     Path to archae_bacteria_plasmid_viral Kraken Database
 -cpu   :i  =>     Number of threads for kraken (def 4)
 -split :i  =>     Split long sequences in this many bp + kmer size (def. 5000).
 -kmer  :i  =>     Kmer used to build database. Fixed at 31.
 -prop  :i  =>     Minimum % proportion of subsequences that needs to be positive (def. 20)
 -stri      =>     Do a stringent search
 -root  :i  =>     TaxID that relates to the 'root' (def 131567)
 -memory    =>     Memory preload the kraken DB (not recommended but sometimes without kraken gets stuck)
 -native    =>     Use kraken-filter. Otherwise use in built filtering

=cut

use strict;
use warnings;
use File::Basename;
use Data::Dumper;
use Pod::Usage;
use Getopt::Long;
use FindBin qw/$RealBin/;
$ENV{PATH} .= ":$RealBin:$RealBin/../kraken/"; 

my $kraken_exec  = &find_program('kraken');
my ($fasta_file,$database_dir,$number_of_seqs,%split_file_data,%average_kmer,$do_stringent,$do_native);
my $kmer_size = 31;
my $memory_map = int(0) ; 
my $cpus = 4;
my $split_size = 5000;
my $prop_cutoff = 20; # perc
my $kmer_filter_cutoff = 0.05;
my $root_taxid = '131567';
GetOptions(
	'infile=s' => \$fasta_file,
	'dir=s' => \$database_dir,
	'cpus:i' => \$cpus,
	'split|size:i' => \$split_size,
	'prop_cut:i' =>\$prop_cutoff,
	'filter_cut:f' => \$kmer_filter_cutoff,
	'memory:i' => \$memory_map,
	'native' => \$do_native,
	'stringent' => \$do_stringent,
	'root:s' => \$root_taxid 
);

pod2usage "No genome file provided" unless $fasta_file && -s $fasta_file;
pod2usage "Invalid Kraken database option\n" unless $database_dir && -d $database_dir && -s "$database_dir/database.kdb" && -s "$database_dir/database.idx";

my $base_out = $fasta_file.".outcomes.".basename($database_dir);

if (-s "$fasta_file.split"){
  &check_split("$fasta_file.split");
}else{
  &split_file($fasta_file);
}

die "No sequences split\n" unless $number_of_seqs && $number_of_seqs > 0;
chomp($number_of_seqs);

	print "Processing $number_of_seqs sequences with $kraken_exec and $database_dir\n";
	my $cmd_options = $memory_map && $memory_map > 0 ? ' --preload' : '';
	$cmd_options .= " --db $database_dir --threads $cpus --fasta-input --output $base_out --only-classified-output $fasta_file.split";
	system($kraken_exec.$cmd_options) unless -s "$base_out";
	if ($do_native){
		system($kraken_exec."-filter --threshold $kmer_filter_cutoff --db $database_dir $base_out > $base_out.$kmer_filter_cutoff");
	}else{
		if ($do_stringent){
			&kraken_filter($base_out,1);
		}else{
			&kraken_filter($base_out);
		}
	}
	system($kraken_exec."-translate  --mpa-format --db $database_dir $base_out > $base_out.labels") unless -s "$base_out.labels";
	system($kraken_exec."-translate  --mpa-format --db $database_dir $base_out.$kmer_filter_cutoff > $base_out.$kmer_filter_cutoff.labels") unless -s "$base_out.$kmer_filter_cutoff.labels";
	die "No kraken output!\n" unless -s "$base_out.labels";
	&process_kraken("$base_out.labels");
	&process_kraken("$base_out.$kmer_filter_cutoff.labels");


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
		my $seq = join('',@lines);
		$seq =~s/\s+//g;
		$id=~s/\s+.+$//g;
		if ($id=~s/\.(\d+)$//){
			my $counter = $1;
			$number_of_seqs++;
			$split_file_data{$id}{$1} = length($seq);
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
		$seq=~s/[^ATGCN]+/N/g;
		# replace strings of N with one N to ensure we get the kmer coverage more accurately
		$seq=~s/N{2,}/N/g;
		my $seq_length = length($seq);
		my @seq_array = split('',$seq);
		my $counter = 1;
		for (my $i=0;$i<$seq_length;$i+=$split_size){
			last if $i >= $seq_length;
			my $k = $i + $split_size - 1 + $kmer_size;
			$k = $seq_length-1 if $k > $seq_length-1;
			my $subseq = join('',@seq_array[$i..$k]);
			next if length($subseq) <= $kmer_size;
			print OUT ">$id.$counter ".($i+1).':'.($k+1)."\n$subseq\n";
			$counter++;
			$number_of_seqs++;
			$split_file_data{$id}{$counter} = length($subseq);
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
	my $file = shift || die;
	my $outfile = $file.".results";
	my %results;
	my %tax_lookup;
	open (IN,$file) || die $!.":$file ";
	while (my $ln = <IN>){
		chomp($ln);
		my @data = split("\t",$ln);
		next unless $data[1] || $data[1] eq 'root';
		my $id = $data[0];
                $id=~s/\s+.+$//g;
		if ($id=~s/\.(\d+)$//){
			my $subseq = $1;
			my $taxstr = $data[1];
			my $taxid = $data[2];
			next if $taxstr eq 'root';
			$tax_lookup{$taxid} = $taxstr if !$tax_lookup{$taxid};
			$results{$id}{$taxid}++;
		}
	}
	open (OUT,">$outfile");
	my %print;
	foreach my $id (keys %results){
		foreach my $taxid (keys %{$results{$id}}){
			my $number_of_subscaffolds = scalar(keys %{$split_file_data{$id}});
			my $prop = sprintf("%.4f",$results{$id}{$taxid}/$number_of_subscaffolds)*100;
			my $taxstr = $tax_lookup{$taxid};
			$print{$id}= $id."\t$taxstr\t$taxid"."\t".$results{$id}{$taxid}."\t$prop\t".$number_of_subscaffolds if $prop > $prop_cutoff;
		}
	}
	print OUT "ID\tTAXONOMY\tTAXID\tSUBSEQ_HIT_THIS\tSUBSEQ_HIT_THIS_%\tTOTAL_SUBSEQS\tAVERAGE_KMERS_HITANYTHING_%\n";

	foreach my $id (sort{$a cmp $b} keys %print){
		my $av = $average_kmer{$id} ? $average_kmer{$id}  : 'N/A';
		print OUT $print{$id}."\t$av\n";
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



sub kraken_filter(){
	# C       flattened_line_2.49     131567  5031    0:3730 2:1 0:1 1:39 131567:1 1434111:1 0:1228
	my $file = shift;
	my $exclude_ambiguous=shift;
	my (%kmer_prop_found,%average_kmer_prop_found);
	open (IN,$file)||die $!;
	open (OUT,">$file.$kmer_filter_cutoff");
	while (my $ln = <IN>){
		next unless $ln && $ln=~/^C\t/;
		chomp($ln);
		my @data = split("\t",$ln);
		my $subscaff_id = $data[1];
		my $taxid = $data[2];
		my $subscaff_length = $data[3];
		my $tax_string = $data[4];
		my @tax_array = split(" ",$tax_string);
		my $kmers_found=int(0);
		my $total_kmers = $subscaff_length - $kmer_size + 1;
		foreach my $tax_str (@tax_array){
			my @tax_data = split(':',$tax_str);
			next unless $tax_data[1];
			if ($exclude_ambiguous){
				$kmers_found += $tax_data[1] if ($tax_data[0] ne '0' && $tax_data[0] ne $root_taxid && $tax_data[0] ne 'A');
			}else{
				$kmers_found += $tax_data[1] if ($tax_data[0] ne '0' && $tax_data[0] ne $root_taxid);
			}
			
		}
		my $prop_kmers = sprintf("%.4f",($kmers_found/$total_kmers));
		$kmer_prop_found{$data[1]}{$data[2]} = $prop_kmers;
		print OUT $ln."\t".($prop_kmers*100)." %\n" if $prop_kmers >= $kmer_filter_cutoff;
	}
	close IN;
	close OUT;

	# too tired to do this properly
	my %counts;
	foreach my $id (keys %kmer_prop_found){
		my $count = int(0);
		foreach my $tax (keys %{$kmer_prop_found{$id}}){
			$average_kmer_prop_found{$id} += $kmer_prop_found{$id}{$tax};
			$count++;
		}
		$average_kmer_prop_found{$id} = ($average_kmer_prop_found{$id} && $count > 0) ? sprintf("%.2f",(($average_kmer_prop_found{$id}/$count) * 100) ) : int(0);
		my $scaffold_id = $id;
		if ($scaffold_id=~s/\.(\d+)$//){
			my $subseq_id = $1;
			$average_kmer{$scaffold_id} += $average_kmer_prop_found{$id};
			$counts{$scaffold_id}++;
		}
	}

	foreach my $scaffold_id (%average_kmer){
		$average_kmer{$scaffold_id} = ($counts{$scaffold_id} && $average_kmer{$scaffold_id} && $average_kmer{$scaffold_id} > 0) ? sprintf("%.2f",($average_kmer{$scaffold_id} / $counts{$scaffold_id})) : 'N/A';

	}

	return "$file.$kmer_filter_cutoff";
}
