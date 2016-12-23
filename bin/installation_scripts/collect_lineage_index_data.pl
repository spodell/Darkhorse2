#!/usr/bin/perl 
# collect_lineage_index_data.pl 
# Sheila Podell 
# June 3, 2007

# this list uses pre-loaded MySQL tables (db_seq_ids and precalc_lineages) from the
# darkhorse database to generate tab-delimited output, containing these fields:

# 	full_id   ncbi_tax_id   species_name   lineage	num_terms

# This version gets relational database connection information from an external
# configuration file (darkhorse.cfg)

# It can select only those lines with lineages having more (good) 
# or less (missing) than a given number of lineage terms  
# the cutoff ($min_lineage_terms) is specified in darkhorse.cfg file

# outputs the list to STDOUT

use warnings;
use strict;
use Getopt::Long;
use DBI;

# get command line arguments
	my ($db, $config_filename, $missing_lineages_only, $good_lineages_only, $all_lineages);
	my $USAGE;
	my $message = qq(
	Usage: $0 -d database_name -c config filename <options>
	
	Options: 
		-m missing lineages only
		-g good lineages only
		-a all lineages (default)
	)."\n\n";
	GetOptions( "d=s" => \$db,
	            "c=s" => \$config_filename,
			    'm' => \$missing_lineages_only,
			    'g' => \$good_lineages_only,
			    'a' => \$all_lineages
				);	     
	
# get SQL database parameters from config file
    my %config_params= &read_config_file($config_filename);
	my $db_name = $config_params{db_name} || "not found";
	my $db_program = $config_params{db_program} || "not found";
	my $db_host = $config_params{db_host} || "not found";
	my $db_user = $config_params{db_user} || "not found";
	my $db_user_password = $config_params{db_user_password} || "not found";
	my $max_lines_per_packet = $config_params{max_lines_per_packet} || "not found";
	my $min_lineage_terms =$config_params{min_lineage_terms} || "not found";	
		
	&check_errs();

# connect to database
	my $db_path = "DBI:$db_program:$db_name:$db_host:";
	my $dbh = DBI->connect($db_path, $db_user, $db_user_password) or die "Cannot connect: $DBI::errstr";

# get all of the species names
	my $sql = qq(SELECT full_id, db_seq_ids.ncbi_tax_id, db_seq_ids.species_name, lineage, num_terms 
	FROM db_seq_ids, precalc_lineages 
	WHERE db_seq_ids.species_name = precalc_lineages.species_name;);

# Prepare and execute the SQL SELECT statement.
    my $sth = $dbh->prepare($sql)
    	or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
	$sth->execute()
    	or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");
	
# Get the lineage information into tab form
	my @row = ();
	my $element = ();
	
	while (@row = $sth->fetchrow_array)
	{
		my $num_lineage_terms = $row[4];			
		if ($num_lineage_terms < $min_lineage_terms) # bad lineages
		{
			if (defined $missing_lineages_only || defined $all_lineages)
			{
				foreach $element (@row)
				{
					print  "$element\t";
				}
				print "\n";
			}	
		}
			
		if ($num_lineage_terms >= $min_lineage_terms)	# good lineages
		{
			unless (defined $missing_lineages_only)
			{
				foreach $element (@row)
				{
					print  "$element\t";
				}				
				print "\n";
			}	
		}			
	}
	
# run a second test to look for entries no species name (primarily from pdb)
# these would not be found in previous join with precalc_lineages;

	$sql = qq(select full_id, db_seq_ids.ncbi_tax_id from db_seq_ids 
	where length(species_name) < 3;);
	
# Prepare and execute the SQL SELECT statement.
    $sth = $dbh->prepare($sql)
    	or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
	$sth->execute()
    	or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");

# write this additional information to ouput	
	while (@row = $sth->fetchrow_array)
	{		
		my $num_lineage_terms = 0;
		if (defined $missing_lineages_only || defined $all_lineages)
		{
			my $count = 0;
			my $lineage = "NULL";
			my $species = "NULL";
			foreach $element (@row)
			{
				$count++;
				print "$element\t";
			}
			
			print "$species\t$lineage\t$num_lineage_terms\n";
		}	
	}
	
# cleanup
	$sth->finish();
	$dbh->disconnect();
	
######################################
# SUBROUTINES
######################################

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

sub check_errs
{	
	if($USAGE || !$db) 
	{
		print STDERR $message;
		exit(0);
	}
	
	if (defined $missing_lineages_only && defined $good_lineages_only)
	{
		print STDERR "Can't have missing_lineages_only AND good_lineages_only\n";
		print STDERR "Make up your mind!\n";
		exit(0);		
	}	
	
	if (defined $all_lineages)
	{
		$good_lineages_only = undef;
		$missing_lineages_only = undef;
		#print STDERR "\nCollecting all lineages\n";
		return;
	}
	elsif (defined $missing_lineages_only)
	{
		$good_lineages_only = undef;
		$all_lineages = undef;
		print STDERR "\nCollecting missing lineages only\n";
		return;
	}	
	elsif (defined $good_lineages_only)
	{
		$missing_lineages_only = undef;
		$all_lineages = undef;
		print STDERR "\nCollecting good lineages only\n";
		return;
	}
	else
	{
		$all_lineages = "yes";
		$good_lineages_only = undef;
		$missing_lineages_only = undef;
		print STDERR "\nCollecting all lineages\n";
		return;
	}
}
__END__

