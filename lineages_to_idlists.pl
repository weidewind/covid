#!/usr/bin/perl
use Getopt::Long;
use IPC::System::Simple qw(capture);
use strict;


my $lineages_file;
my $pangolin_file;
my $output;
my $log;

GetOptions (	
		'lineages_file=s' => \$lineages_file,
		'pangolin_file=s' => \$pangolin_file,
		'output=s' =>\$output,
		'log=s' =>\$log,
		);
$| = 1;
open P, "<$pangolin_file" or die "Cannot open $pangolin_file : $! \n";
my @pango = <P>;
close P;

open L, "<$lineages_file" or die "Cannot open $lineages_file : $! \n";
my @lins = <L>;
close L;

open LOG, ">>$log" or die "Cannot open $log for writing : $! \n";
print LOG "Number of ids corresponding to clade, according to pangoline.all\n";

@lins = sort { scalar split(/\./,$b) <=> scalar split(/\./, $a) } @lins;
foreach my $l(@lins){
	$l =~ s/[\n\r]+//;
	my $lname = $l;
	print($l."\n");
	my @idlist;
	if (substr($l,length($l)-1,length($l)) eq '$'){
		$l = substr($l, 0, length($l)-2);
		@idlist = grep {/\t\Q${l}\E$/} @pango; #take ${l} literally and not as regex
		@pango = grep {!/\t\Q${l}\E$/} @pango; 
	}
	else {
		@idlist = grep {/\t\Q${l}\E$/ || /\t\Q${l}\E\./} @pango;
		@pango = grep {!(/\t\Q${l}\E$/ || /\t\Q${l}\E\./)} @pango; 
	}
	my $outpath = $output.$lname;
	open OUT, ">$outpath" or die "Cannot open $outpath: $!\n";
	print LOG $lname."\t".scalar @idlist."\n";
	foreach my $str(@idlist){
	#	my @sp = split(/\t/,$str);
	#	print OUT $sp[0]."\n";
		print OUT $str."\n";
	}
	close OUT;
}
close LOG;