#!/usr/bin/perl 
# cluster_lineages.pl
# Sheila Podell
# February 18, 2017

# takes as input a file of lineages (e.g. from a metagenome)
# counts abundance of lineage strings
# filters out those below given minimum percentage of whole (default = 5%)
# collapses shorter lineages with same number of hits into longer ones

# outputs tab-delimited list (ordered by descending number of members)
	# nickname	num_members	percentage	num_terms	cluster
	
use warnings; 
use strict;
use Getopt::Long;

# Global defaults
	my ($dh_smry_filename, $ctg_index_filename);
	my $min_pct = .05;
	
# get command line arguments
	my $USAGE;
	my $message = qq(
	Usage: $0 -i darkhorse_smry_file 
	optional parameter: -n min_pct_matches \(default = 0.05 (5%)\)	
	)."\n";
	
	GetOptions( "i=s" => \$dh_smry_filename,
"n=f" => \$min_pct,	
        );
    
# validate input
	if($USAGE || !$dh_smry_filename) 
    {
		print STDERR "$message";
		exit(0);
	}
	unless (defined $min_pct 
		&& $min_pct =~ /[\d\.].+$/
		&& $min_pct >0
		&& $min_pct < 1 )
	{
		print STDERR "   ERROR: min_pct_matches must be a positive value between 0-1.";
		exit(0);
	}
	
#  get clusters
	my %lineages = ();	# key = lineage name # value = object reference
	my $current_lineage_object;
	my %clusters = (); # key =  cluster name # value = object reference
	my $current_cluster;
	my $input_lines = 0; 
	open (DH_SMRY,"$dh_smry_filename") or die "can't open input file $dh_smry_filename\n$!\n";		
	while (<DH_SMRY>) 
	{
		next if ($_ =~ /^\s/);
		next if ($_ =~ /query/); #skip header line
		chomp;
		$input_lines++;
		my @tmp = split "\t", $_;
		my $lineage_string = $tmp[14];		
		unless (defined $lineage_string && exists $lineages{$lineage_string})
		{ 
			$current_lineage_object = &new_lineage("lineage", "$lineage_string");
			$lineages{$lineage_string} = $current_lineage_object;		
		}
	
	# populate clusters
		my @cluster_terms = split ";" , "$lineage_string";
		my $cluster_name = "$cluster_terms[0]";
		for (my $i = 0; $i < scalar @cluster_terms; $i++)
		{					
			my $next_term = $cluster_terms[$i];
			# create cluster name by joining all terms up to and including current term
			my $previous_term = $cluster_terms[$i-1];
			if (scalar @cluster_terms > 1)
			{				
				$cluster_name = join ";", (@cluster_terms[0..$i]);
			} 
					
			unless (exists $clusters{$cluster_name})
			{ 
				$current_cluster = &new_cluster("cluster", $cluster_name);
				$clusters{$cluster_name} = $current_cluster;
			}
			
			$clusters{$cluster_name}->{num_members}++;
		}				
	}
	close DH_SMRY;	
	
# create list of clusters ordered by number of matches
	my @sorted_cluster_objects = sort {
		$b->{num_members} <=> $a->{num_members}
		||
		$a->{num_terms} <=> $b->{num_terms}		
	} values %clusters;
	
# check to see any clusters are a subset with exactly the same number of matches
	my @sorted_clusters = ();
	my $current_cluster_object = "";
	my $current_name;
	my $current_num_matches = 0;
	my $prev_cluster_object = "";
	my $prev_name = "";
	my $prev_num_matches = 0;

	for (my $i = 0; $i < scalar @sorted_cluster_objects; $i++)
	{
 		$current_cluster_object = $sorted_cluster_objects[$i];
 		$current_name = $current_cluster_object->{name};
 		$current_num_matches = $current_cluster_object->{num_members};
 		#print STDERR "current name=$current_name print STDERR   current_num_matches=$current_num_matches\n";
		
		if ($i > 0)
		{
			$prev_cluster_object = $sorted_cluster_objects[$i-1];
			$prev_name = $prev_cluster_object->{name};
			$prev_num_matches = $prev_cluster_object->{num_members};			
		}
		else
		{
			$prev_name = "blank";
			$prev_num_matches = 0;
		}
				
	# get rid of previous object if has same # matches and is subset of current	
		if ( $i >0
			&& $current_num_matches == $prev_num_matches 
			&& $current_name =~ /^$prev_name.+/ )
		{			
			pop @sorted_clusters;		
		}		

		push @sorted_clusters, $current_cluster_object->{name};
	} 

# output results
	my @headers = (
	"nickname",
	"num_members",
	"percentage",
	"num_terms",
	"cluster",
	);
	my $header_string =  join "\t", @headers;
	print "$header_string\n";
	
	foreach my $clstr_name (@sorted_clusters)
	{
		my $num_members = $clusters{$clstr_name}->{num_members};
		my $nickname = 	$clusters{$clstr_name}->{nickname};
		my $num_terms =  $clusters{$clstr_name}->{num_terms};
		
		my $pct = $num_members/$input_lines;
		next unless $pct > $min_pct;
	
	# consolidate
		
		my $formatted_pct = 100*$pct;
		$formatted_pct = sprintf("%0.0f", $formatted_pct);
		$formatted_pct = "$formatted_pct"."%";
		
		my @print_cols = ($nickname, $num_members, $formatted_pct, $num_terms,  $clstr_name, );
		my $printstring = join "\t", @print_cols;
		print "$printstring\n";
	}

#####################
# SUBROUTINES
#####################

sub new_lineage {
  my ($className, $lineage_name, $param) = @_;
  my $self = {};
  bless $self, $className;
	$self->{name} = $lineage_name;
	$self->{num_matches} = 0;
	my @tmp = split ";" , $lineage_name;
	$self->{num_terms} = scalar @tmp;

  return($self);
}

sub new_cluster {
  my ($className, $cluster_name, $param) = @_;
  my $self = {};
  bless $self, $className;
	$self->{name} = $cluster_name;
	my @tmp = split ";", "$cluster_name";
	$self->{nickname} = $tmp[-1]; # last term
	$self->{num_members} = 0;
	$self->{num_terms} = scalar @tmp;

  return($self);
}


      
