#!/usr/bin/perl
# count_children.pl
# Sheila Podell
# February 2, 2017

# Using DarkHorse2 database, measures taxonomic diversity of database representatives 
# for a given taxonomic group 

# does this by counting distinct terms immediately adjacent (to the right) in lineage strings e.g.
	# Bacteria;Proteobacteria;Deltaproteobacteria;Syntrophobacterales;Syntrophaceae;Syntrophus
	# Bacteria;Proteobacteria;Deltaproteobacteria;Myxococcales;Sorangiineae;Polyangiaceae;Chondromyces
	# Bacteria;Planctomycetes;Planctomycetia;Planctomycetales;Planctomycetaceae
 
# input: file containing list of query terms (one per line) from a lineage string
# example (don't use # sign before!):
# 	Archaea
#	Bacteria
#       Proteobacteria
#	Deltaproteobacteria
	
# output: tab-delimited count for number of distinct terms found one position to the right
# (direct descendent children of the test term, but not grandchildren, etc) 
#       Archaea	    0
# 	Bacteria	2
#	Proteobacteria	1
#	Deltaproteobacteria	2

use warnings;
use strict;
use Getopt::Long;
use DBI; 

# get command line arguments
	my ($db, $taxonomy_term_listfile, $lineage_listfile, $sth, $sql, $config_filename, $printline, $debug);
	my $USAGE;
	my $message = qq( 
	Usage: $0  -c config_filename -t query_listfile <optional: -l lineage_listfile)."\n\n";
	GetOptions( 
			    "t=s" => \$taxonomy_term_listfile,
			    "c=s" => \$config_filename,
			    "l=s" => \$lineage_listfile,
			    "d=s" => \$debug
				);	
				
# check arguments
	if($USAGE || !$config_filename  || !$taxonomy_term_listfile)
	{
		print STDERR $message;
		exit(0);
	}
	unless (-s $taxonomy_term_listfile)
	{
		print STDERR "  Couldn't find file containing taxonomy term listfile $taxonomy_term_listfile \n";
		exit(0);
	}	
	unless (-s $config_filename)
	{
		print STDERR "  Couldn't find config file $config_filename \n";
		exit(0);
	}	
	
# get SQL database parameters from config file
	my %config_params = &read_config_file($config_filename);
	my $db_name = $config_params{db_name} || "not found";
	my $db_program = $config_params{db_program} || "not found";
	my $db_host = $config_params{db_host} || "not found";
	my $db_user = $config_params{db_user} || "not found";
	my $db_user_password = $config_params{db_user_password} || "not found";
	my $max_lines_per_packet = $config_params{max_lines_per_packet} || "not found";
	my $min_lineage_terms = $config_params{min_lineage_terms} || "not found";
	
# connect to database
	my $db_path = "DBI:$db_program:$db_name:$db_host:";
	my $dbh = DBI->connect($db_path, $db_user, $db_user_password) or die "Cannot connect: $DBI::errstr";
	my @row = ();
	my $element;

# get list of input names, put into hash
	open (INPUT, $taxonomy_term_listfile)  or die "can't input file $taxonomy_term_listfile\n$!\n";
	my %input_tally = (); #key = term value = num children
	my @input_order =();  # keep track so output can be in same order as inpt
	while (<INPUT>)
	{  	
		$_ =~ s/\r\n/\n/s; # convert windows to unix line endings, if necessary
		next if $_ =~ /^\s*$/;	# ignore blank lines
		chomp;
		if (defined $_ && length $_ > 1)
		{		
			unless (exists $input_tally{$_}) # don't use  duplicates
			{
				push @input_order, $_;
				$input_tally{$_} = 0;
			}
		}    
	}    
	close INPUT;
	
	my $num_input_tally = scalar (keys %input_tally);
  	if ($num_input_tally  > 0)
  	{
  		print STDERR "  Found $num_input_tally taxonomy terms to test\n";
  	}	
	else
	{
		print STDERR "  ERROR: no taxonomy query terms founn.";
		exit(0);
	}
	
# if no listfile provided, generate one	
	my %lineage_list = ();
	my %uniq_children = ();
	
	unless (defined $lineage_listfile && -s $lineage_listfile)
	{
		print STDERR "  Generating lineage list by SQL query from $db_name\n";
		print STDERR "  This may take a while . . .\n";
		
	#perform sql query to get all lineage candiates, put into listfile
		$lineage_listfile = "$db_name"."_uniq_lineages";
		open (OUTPUT, ">$lineage_listfile") or die "couldn't open lineage listfile $lineage_listfile for writing, $!\n";
		$sql  = qq (
select distinct lineage
from lineage_index 
where num_terms > 1;		
);
	$sth = $dbh->prepare($sql)
		or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
	$sth->execute()
	or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");
	
	my @row = ();
	my $element = ();
	
	# write immediately to file - no buffer wait (perl cookbook 7.12, p. 248) 
	# that way, in case database disconnects early, you've still printed results
	# up to that point.
	
	$| = 1;
			
	while (@row = $sth->fetchrow_array)
	{									
		my ($lineage_name) = $row[0];
		chomp $lineage_name;
		print OUTPUT "$lineage_name\n";			
	}
	$dbh->disconnect();
	close OUTPUT;
}

# parse the lineage listfile	
	if (defined $lineage_listfile && -s $lineage_listfile)
	{
		print STDERR "  Using lineage listfile $lineage_listfile\n";
		open (INPUT, $lineage_listfile)  or die "can't open lineage listfile $lineage_listfile\n$!\n";
		while (<INPUT>)
		{
			$_ =~ s/\r\n/\n/s; # convert windows to unix line endings, if necessary
			next if $_ =~ /^\s*$/;	# ignore blank lines
			next if ($_ =~ /^#/);  #ignore commented lines
			chomp;
			if (defined $_ && length $_ > 3)
			{
				next unless $_ =~ /;/;  # make sure there are at least 2 terms
				$lineage_list{$_}++;
			} 
			my @lineage_terms = split ";", $_;  	
			# go through each of tax terms, count child hits;
			foreach my $input_query (@input_order)
			{
				#skip lines without term 
				next unless $_ =~ /$input_query/;
				# find position in lineage array
				my $lineage_position = -1;
				foreach my $next (@lineage_terms)
				{
					$lineage_position++;
					if ($next =~ /$input_query/)
					{
						my $child_position = $lineage_position + 1;
						my $child_name = $lineage_terms[$child_position];
						next unless (defined $child_name);
						$uniq_children{$child_name}++;
						unless ($uniq_children{$child_name} >1)
						{
								$input_tally{$input_query}++;
						}						
					}
				}				
			}
		}
		close INPUT;
	}
	
	my $num_distinct_lineages = scalar keys %lineage_list;
	
	if ($num_distinct_lineages > 0)
	{
		print STDERR "  Searched $num_distinct_lineages distinct lineages for child names.\n";
	}
	
	
	
	
	# print output
	foreach my $term (@input_order)
	{
		print "$term"."\t"."$input_tally{$term}\n"; ;
	}	
	

	
################################################
# SUBROUTINES
################################################
sub read_config_file
{
    my ($filename) = @_;
    open (INPUT, $filename)  or die "can't open DarkHorse config file $filename\n$!\n";
	my %dh_config = ();
	while (<INPUT>)
	{
		next if ($_ =~ /^\s+$/);
		next if ($_ =~ /^#/);
		if ($_ =~ /\[(.+)\]=\s*(.+)/)
		{
			$dh_config{$1}=$2;		
		}			
	}
	close INPUT;			
    return %dh_config;
}


	
	
