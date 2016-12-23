#!/usr/bin/perl 
# merge_annot_protid_taxid_tables.pl 
# Sheila Podell 
# November 14, 2016

# this list uses pre-loaded MySQL tables (annotation and protid_taxid) from the
# darkhorse database to generate tab-delimited output, containing these fields:

# 	full_id   ncbi_tax_id   species_name   seq_length	annotation

# This version gets relational database connection information from an external
# configuration file (darkhorse.cfg)

# outputs the list to STDOUT

use warnings;
use strict;
use Getopt::Long;
use DBI;

# get command line arguments
	my ($db, $config_filename, $missing_lineages_only, $good_lineages_only, $all_lineages);
	my $USAGE;
	my $message = qq(
	Usage: $0  -c config filename 
	
	)."\n\n";
	GetOptions( 
	            "c=s" => \$config_filename,
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
	
	my $sql = qq(select protid_taxid.full_id, ncbi_tax_id, species_name, annotation, seq_len
	from protid_taxid, annotation
	where protid_taxid.full_id = annotation.full_id;);

# Prepare and execute the SQL SELECT statement.
    my $sth = $dbh->prepare($sql)
    	or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
	$sth->execute()
    	or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");

# write this additional information to ouput	
	my @row = ();
	while (@row = $sth->fetchrow_array)
	{		
		my $row_txt = join "\t", @row;		
		print "$row_txt\n";
		@row = ();
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
	if($USAGE) 
	{
		print STDERR $message;
		exit(0);
	}
}
__END__

