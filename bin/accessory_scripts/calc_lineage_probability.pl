#!/usr/bin/perl
# calc_lineage_probability.pl
# Sheila Podell 
# June 2, 2007

# takes as input a file containing  rows of colon delimited taxonometric lineages
# from ncbi. example: 
#	Bacteria; Cyanobacteria; Nostocales; Nostocaceae; Nostoc.

# Gets a hierarchy number for each term (is it at top or bottom of lineage?) 
# counts number of times each term is present in the genome

# sends tab delimited result to STDOUT. Fields:
# 	lineage 	LPI

# revised October 1, 2013 to accomodate lineage strings where the same term occurrs
# more than once (e.g. Actinobacteria, which is both a phylum and a class name)

# revised November 5 2016 to output lineage identical to input (restore internal spaces)

use warnings; 
use strict;

#	Get list of lineages from file specified on command line
	unless (@ARGV == 1)
	{
		print "Usage: $0 list_file\n";
		exit (0);
	}
	my $infilename = $ARGV[0];
	open (INPUT, "$infilename") or die "can't open input file $infilename \n$!\n";
	
# get lineage terms into data structures to tally frequencies for whole data set
	my ($current, $name) = ""; 		# object reference 		
	my %terms = ();			# key = term name  value = object reference
	my $rowcount = 0;
	while (my $line = <INPUT>) 
	{
		next if $line =~ /^\s+$/;	#skip blank lines
		chomp $line;
		
		$line =~ s/\.//g;	 		#get rid of period at end
		$line =~ s/ //g;			#get rid of extra spaces
		my @tmp = split /\;/, $line;
		for (my $i = 0; $i < scalar @tmp; $i++)
		{			
			$name = $tmp[$i];													
			unless (exists $terms{$name})
			{
				my @parameters = ($name, $i, $rowcount);
				$current = &new_term("term", \@parameters);
				$terms{$name} = $current;
			}
			$terms{$name}->{count}++;
			$terms{$name}->{hierarchy_sum} += $i;
		}
		$rowcount++;
	}
	close INPUT;

# figure out average hierarchical position of each term
# for giving higher weights to more general taxonomic terms (e.g. kingdom versus species)
	my ($num_hits, $sum, $avg_hierarchy);
	foreach (keys %terms)
	{
		$num_hits = $terms{$_}->{count};
		$sum = $terms{$_}->{hierarchy_sum};
		$avg_hierarchy = int($sum /$num_hits);
		$terms{$_}->{avg_hierarchy} = $avg_hierarchy;
	}

# find out how many entries there are for each hierarchical category
	my %categories = ();
	my $level = "";
	
	foreach (keys %terms)
	{
		$level =  $terms{$_}->{avg_hierarchy};
		$categories{"$level"} += 	$terms{$_}->{count};	
	}	

# sorted by hierarchical category, decreasing frequency
	foreach (sort 
		{
			$terms{$a}->{avg_hierarchy} <=> $terms{$b}->{avg_hierarchy}
			||
			$terms{$b}->{count} <=> $terms{$a}->{count}
		} 
		keys %terms)
	{
		
	# calculate probability for term within it's hierarchical category
		my $frequency = $terms{$_}->{count};
		my $current_category = $terms{$_}->{avg_hierarchy};
		my $level_count = $categories{$current_category};
		$terms{$_}->{probability} = $frequency/$level_count;
	# calculate weight - add 1 to current category to avoid dividing by zero	
		$terms{$_}->{weighted_probability} = ($frequency/$level_count)/($current_category+1);
	}	
		
# append probabilities to original input file
	open (INPUT, "$infilename") or die "can't open input file $infilename \n$!\n";
	
# calculate composite probabilty for each lineage
	while (my $line = <INPUT>) 
	{
		next if $line =~ /^\s+$/;	#skip blank lines
		chomp $line;
		my $input_lineage = $line;
		$line =~ s/\.//g;	 		#get rid of period at end
		$line =~ s/ //g;			#get rid of extra spaces
		
		my @tmp = split /\;/, $line;
			
		my $add_string = "";
		my $num_terms = scalar @tmp;		
		my $sum = 0;
		my $max_possible = 0;
		my %current_terms = ();

		for (my $i = 0; $i < scalar @tmp; $i++)
		{						
			my $next_term = $terms{$tmp[$i]};
		# deal with special case where same name used for both Phylum and Class (e.g. Actinobacteria)
			$current_terms{$next_term}++;
			next if ($current_terms{$next_term} > 1);
						
			$add_string .= "\t$next_term->{weighted_probability}";
			$sum += $next_term->{weighted_probability};
			$max_possible += 1/($i+1);
		}
# normalize sum to a probability - divide by maximum sum possible for # terms
# make sure that output lineage matches initial input
		my $probability = $sum/$max_possible;	
		print "$input_lineage$probability\n";		
	}
	close INPUT;

print "\n";

############################################
# SUBROUTINES
############################################	

sub new_term {
  my ($className, $param) = @_;
  my $self = {};
  bless $self, $className;
  my @properties = @$param;
  $self->{name} = $properties[0];  
  $self->{hierarchy_sum} = $properties[1];
  $self->{input_order} = $properties[2];
  $self->{count} = 1;
  $self->{avg_hierarchy} = "";
  $self->{probability} = "";

  return($self)
}

sub term_toString {

	my ($term) = @_;
	my $string = qq(
$term
name = $term->{name}
hierarchy_sum = $term->{hierarchy_sum}; 
input_order = $term->{input_order}
count = $term->{count}
avg_hierarchy = $term->{avg_hierarchy}

);
print STDERR $string;

}
