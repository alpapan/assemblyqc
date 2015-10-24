#!/usr/bin/env perl


use strict;
use warnings;
use File::Basename;

my $auto;


if ($auto){
	my (%list);
	my @get = ("ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/plasmid.1.1.genomic.fna.gz","ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/plasmid.2.1.genomic.fna.gz");
	$list{'bacteria'} = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt';
	$list{'viral'} = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt';
	$list{'archae'} = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt';
	foreach my $tax (keys %list){
		mkdir($tax) unless -d $tax;
		chdir($tax);
		system("wget -c -N ".$list{$tax});
		&process("assembly_summary.txt");
		chdir("..");
	}
	foreach my $file (@get){
		my $existing = $file;
		$existing=~s/\.gz$//;
		if (!-s $existing){
			system("wget -c -N $file");
			system("gunzip $file");
		}
	}
}

else{

my $txt = shift;
my $do_contig = shift;
 &process($txt,$do_contig);

}


###############################
sub process(){
my ($txt,$do_contig) = @_;
die unless $txt && -s $txt;

print "Processing $txt\n";
my %hash;
open (TXT,$txt)||die;
my $header = <TXT>;
while (my $ln=<TXT>){
	chomp($ln);
	my @data = split("\t",$ln);
	next unless $data[19];
	my $file = basename($data[19])."_genomic.fna";
	next unless $data[10] eq 'latest';
	next unless $data[11] eq 'Complete Genome' || $data[11] eq 'Scaffold' || $data[11] eq 'Chromosome' || ($data[11] eq 'Contig' && $do_contig);
	if (!-s $file && !-s $file.".gz"){
		system("wget -c -N ".$data[19]."_genomic.fna.gz");
	}

	$hash{$file}=$data[5];
}
close TXT;

my @fasta = glob("*fna");
push(@fasta,glob("*fna.gz"));
open (OUT,">all_seqs.tax");
my $orig_sep = $/;
$/ = ">";
my $counter;
my $total = scalar(@fasta);
foreach my $fsa (@fasta){
	if ($fsa=~s/\.gz$//){
		system("gunzip $fsa.gz");
	}
	$fsa = basename($fsa);
	open (FSA,$fsa)||die;
	while (my $record=<FSA>){
		chomp ($record);
		next unless $record;
		my @lines = split("\n",$record);
		my $id = shift @lines;
		if ($id=~/^(\S+)/){
			my $taxid = $hash{$fsa};
			unless ($taxid){
				warn "\nNo tax for $fsa\n";
				next;
			}
			print OUT ">$1|kraken:taxid|$taxid\n".join("\n",@lines)."\n";
		}
	}
	close FSA;
	$counter++;
	print "Processed $counter / $total    \r";
}
$/ = $orig_sep;
close OUT;
}
