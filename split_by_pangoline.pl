#!/usr/bin/perl
use Getopt::Long;
use IPC::System::Simple qw(capture);
use strict;


my $fasta;
my $meta;
my $output;

GetOptions (	
		'fasta=s' => \$fasta,
		'meta=s' => \$meta,
		'output=s' =>\$output,
		);


		
my %filehandles;
my %counter;
my $total;
open META, "<$meta" or die "Cannot open $meta: $!";
while (<META>){
	my @splitter = split(/\t/); #strain  virus   gisaid_epi_isl  genbank_accession       date    region  country division     location        region_exposure country_exposure        division_exposure   segment  length  host    age     sex     Nextstrain_clade        pangolin_lineage    GISAID_clade     originating_lab submitting_lab  authors url     title   paper_url   date_submitted   purpose_of_sequencing
	my $id = $splitter[2];
	my $pango = $splitter[18];
	unless (exists $counter{$pango}) {$counter{$pango} = 0;}
	$counter{$pango} += 1; 
	$total += 1; 
	# if (exists $filehandles{$pango}){
		# my $filehandle = $filehandles{$pango};	
		# print $filehandle $id;
	# }
	# else {
		# print "Found pangoline $pango\n";
		# local *FILE;
		# my $path  = $output."_".$pango."_ids";
		# open FILE, ">$path" or die "Cannot create $path";
		# FILE->autoflush(1);
		# print FILE $id;
		# $filehandles{$pango} = *FILE;
	# }
}
close META;

my $log = $output.".log";
open LOG, ">$log" or die "Cannot open $log for writing: $!\n";
print LOG "Total $total entries in $meta\n";
my $facount = capture ('grep -c ">" '.$fasta);
print LOG "Total $facount entries in $fasta\n";
foreach my $pango (sort { $counter{$a} <=> $counter{$b} } keys %counter){
	print LOG $pango." ".$counter{$pango}."\n";
}
close LOG;

# foreach my $pango (keys %filehandles){
	# print "Creating fasta from pangoline $pango ..\n";
	# close $filehandles{$pango};
	# my $ids  = $output."_".$pango."_ids";
	# my $out = $output."_".$pango.".fasta";
	# my $logs = capture ('faFilter -namePatList='.$ids.' '.$fasta.' '.$out);
	## to do: inverse filter, id not in composite list: -v 
 	# print $logs."\n";
# }