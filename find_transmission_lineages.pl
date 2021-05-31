#!/usr/bin/perl

# Input: 1) --tree newick tree 2) --states fasta with states for all nodes (including leaves)
# Output format: list of nodes closest to root where probability of Russian state (2) turns to 1.0
# >Entry_node
# list of strains in the corresponding subtree (if entry_node is not a leaf)

use Getopt::Long;
use Bio::Phylo::IO;
use strict;

my $tree;
my $states;
my $output;

GetOptions (	
		'tree=s' => \$tree,
		'states=s' => \$states,
		'output=s' =>\$output,
		);
$| = 1;

print "Parsing tree.. \n";
my $st = time();
my $tree = parse_tree($tree);
my $end = time();
my $runtime = $end-$st;
print "Finished in $runtime sec \n";
my %states;
print "Parsing states.. \n";
$st = time();
open ST, "<$states" or die "cannot open states file $states $! \n";
while(<ST>){
	my ($nname, $state1_prob, $state2_prob) = split(/\s+/);
	$states{$nname} = $state2_prob; # probability of being russian
	print "nname $nname state $state2_prob\n";
}

close ST;
$end = time();
$runtime = $end-$st;
print "Finished in $runtime sec \n";

my %lineages; # lineages{entry_node_name} = @strains

print "Looking for lineages.. \n";
$st = time();
my @array;
my_visit_depth_first ($tree->get_root, \@array,  \&collect_lineage_members,\&not_exported, "");
$end = time();
$runtime = $end-$st;
print "Finished in $runtime sec \n";

foreach my $entry_nname(keys %lineages){
	print ">".$entry_nname."\n";
	print join("\n", @{$lineages{$entry_nname}})."\n";
}

sub collect_lineage_members {
	my $node = $_[0];
	my $entry_nname = $_[1];
	#print "Node ".$node->get_name()." Entry name $entry_nname\n";
	if ($entry_nname){
		if ($node->is_terminal){push $lineages{$entry_nname}, $node->get_name();}
	}
	else {
		#print "Probability of being russian is ".$states{$node->get_name()}."\n";
		if ($states{$node->get_name()} == 1 ){ # probability of being russian is 1
			$entry_nname = $node->get_name();
			print "Lineage started at $entry_nname\n";
			my @linlist;
			$lineages{$entry_nname} = \@linlist;
			if ($node->is_terminal){push $lineages{$entry_nname}, $node->get_name();}
		}
	}
	return $entry_nname;
} 

sub not_exported {
	my $node = $_[0];
	my $entry_nname = $_[1];
	if ($entry_nname && $states{$node->get_name()} == 0){return 0;} # already been to russia but now probability of being russian is 0
	else {return 1;}
}


sub parse_tree {
					my $tree_file = $_[0];
					open TREE, "<$tree_file" or die "Cannot open file ".$tree_file."\n";
					# get a newick string from some source
 					my $tree_string;
 					my $t;
 					while($t = <TREE>){
 						$t =~ s/\n$//;
 						$tree_string = $tree_string.$t;
 					}
 					 close TREE;

 					# Call class method parse from Bio::Phylo::IO
 					# note: newick parser returns 'Bio::Phylo::Forest'
                    # Call ->first to retrieve the first tree of the forest.
 					my $tree = Bio::Phylo::IO->parse(
  					  -string => $tree_string,
   					  -format => 'newick'
 					)->first;

 					return $tree;
	
} 

 sub my_visit_depth_first {
		my $node = $_[0];
		my @array = @{$_[1]};
		my $action_callback = $_[2];
		my $check_callback = $_[3];
		my $action_result = $_[4];
		push @array, $node;
		my $ok = &$check_callback($node, $action_result);
		#print "Processing node ".$node->get_name()."\n";
		#print "Is it ok? (not exported) $ok \n";
		if ($ok) {$action_result = &$action_callback($node, $action_result)}; # current $entry_nname
		if (! $node->is_terminal && $ok){
			my $i = 0;
			while($node->get_child($i)){
				@array = my_visit_depth_first($node->get_child($i), \@array, \&$action_callback, \&$check_callback, $action_result);
				$i++;
			}
		}	
		$node = pop @array;
		#print $node->get_name()."\n";
		return @array;
    }