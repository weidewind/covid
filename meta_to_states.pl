#!/usr/bin/perl

# finds all strains in the tree file and looks them up in meta file
# output: tab-delimemted file, strain\tstate, where state is 1 or 2: 1 - the strain is foreign, 2 - the strain is Russian 

use Getopt::Long;
use IPC::System::Simple qw(capture);
use Spreadsheet::XLSX;
use strict;

my $tree;
my $meta;
my $output;

GetOptions (	
		'tree=s' => \$tree,
		'meta=s' => \$meta,
		'rpn_meta=s' => \$rpn_meta,
		'output=s' =>\$output,
		);
$| =1;
open TREE, "<$tree" or die "Cannot open $tree $! \n";
my @ids;
while(<TREE>){
	my @matches = $_ =~ /EPI_ISL_[0-9]+/g;
	push @ids, @matches;
}
close TREE;

open META, "<$meta" or die "Cannot open $meta $! \n";
my %id_to_country;
while(<META>){
	my @splitter = split(/\t/);
	$id_to_country{$splitter[2]} = $splitter[6];
}
close META;

my $excel = Spreadsheet::XLSX->new($rpn_meta);
my $sheet = $excel -> {Worksheet}-> [0];

open OUT, ">$output" or die "cannot open $output for writing $! \n";
foreach my $id (@ids){
	my $state = $id_to_country{$id} eq "Russia" ? 2 : 1;
	print OUT $id."\t".$state."\t".$id_to_country{$id}."\n";
}
close OUT;