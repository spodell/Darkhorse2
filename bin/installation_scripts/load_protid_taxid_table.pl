#!/usr/bin/env perl
# load_protid_taxid_table.pl
# Sheila Podell
# November 5, 2016

use warnings;
use strict;
use Getopt::Long;
use DBI;

# creates, populates a Mysql table correlating sequence id numbers with ncbi taxonomy id 
# expects three files as input
	# 1) darkhorse_config file
	# 2) list of valid ncbi tax ids - previously determined to meet the following rules
		# present in DarkHorse mysql table "precalc_lineages"
		# tax_id != 0
		# num lineage terms > 2
		# valid species name (not multispecies identified only by genus) 
	
	# 3) tab-delimited input file with 2 fields
		#prot_id    ncbi_tax_id 
	
		# example: ncbi file prot.accession2taxid to middle two columns, e.g.
		#	ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid
		#	LC_ALL=C sort prot.accession2taxid |cut -f2,3 > protid_taxid_input.tab
			#prot_id    ncbi_tax_id
	
use warnings; 
use strict;
use Getopt::Long;
use DBI;

# get command line arguments, defaults 
	my ($datafile, $config_filename, $valid_idlist_filename);
	my ($db, $sth, $sql, $overwrite, $update, $result);
	my $USAGE;
	my $message = qq(
	Usage: $0  -c config filename -i protid_taxid_input_file -v valid_id_list_file)."\n\n";
	GetOptions( 
			    "c=s" => \$config_filename,
			    "i=s" => \$datafile,
			    "v=s" => \$valid_idlist_filename,
			    
			    'o'   => \$overwrite,		#allow table overwrites
			    'u'   => \$update,
				);	     
	
	if($USAGE || !$config_filename || !$valid_idlist_filename) 
	{
		print STDERR $message;
		exit(0);
	} 
	
	unless (-s $config_filename)
	{
		print STDERR "  Couldn't find config file $config_filename \n";
		exit(0);
	}
	
	unless (-s $datafile)
	{
		print STDERR "  Couldn't find data input file $datafile \n";
		exit(0);
	}
	unless (-s $valid_idlist_filename)
	{
		print STDERR "  Couldn't find data input file $valid_idlist_filename \n";
		exit(0);
	}
	
# load valid_idlist	into searchable hash
	my %valid_ids = ();
	open (INPUT,"$valid_idlist_filename") or die "can't open input file $valid_idlist_filename\n $!\n";
	while (<INPUT>)
	{
		chomp;
		next if ($_ =~ /^\s+$/);
		my @tmp = split "\t", $_;
		my $taxid = $tmp[0];
		my $species_name = $tmp[1];
		$valid_ids{$taxid}=$species_name;
	}
	close INPUT;
	
# get SQL database parameters from config file
    my %config_params= &read_config_file($config_filename);
	my $db_name = $config_params{db_name} || "not found";
	my $db_program = $config_params{db_program} || "not found";
	my $db_host = $config_params{db_host} || "not found";
	my $db_user = $config_params{db_user} || "not found";
	my $db_user_password = $config_params{db_user_password} || "not found";
	my $max_lines_per_packet = $config_params{max_lines_per_packet} || "not found";
	
# connect to database
	my $db_path = "DBI:$db_program:$db_name:$db_host:";
	my $dbh = DBI->connect($db_path, $db_user, $db_user_password) or die "Cannot connect: $DBI::errstr"; 	

# check to see if table already exists
	my $tablename = "protid_taxid";
	my $found_table = 0;
	$sql = "show tables;";
	$sth = $dbh->prepare($sql)
    	or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
	$sth->execute()
    	or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");
	my @row = ();
	my $element = ();
	while (@row = $sth->fetchrow_array)
	{
		foreach $element (@row)
		{
			if ($element =~/$tablename/)
			{
				$found_table = 1;
				unless (defined $overwrite || defined $update)
				{
					print STDERR "  Table $tablename already exists.\n";
					exit(0);
				}
				last;
			}
		}
	}
	
# if over-write option specified, drop table from database.
	if (defined $overwrite)
	{
		print STDERR "Over-writing table $tablename.\n";
		$sql = qq(DROP TABLE IF EXISTS `$tablename`);
		$result = $dbh->do($sql);		
	}

	if (defined $update)
	{
		print STDERR "updating table $tablename.\n";
	}

# Compose creation command
	if (defined $overwrite || $found_table == 0)
	{
		$sql = qq(CREATE TABLE `$tablename` (
	   `full_id` varchar(50)  NOT NULL,
	   `ncbi_tax_id` int(11)  NULL,
	   `species_name` varchar(100)  NULL,
		PRIMARY KEY  (`full_id`),
		INDEX (ncbi_tax_id)
	 ) ;\n);
 		
 		$result = $dbh->do($sql); 	
		if ($result)
		{
			print STDERR "New table created.\n";
		}
		else
		{
			print STDERR "Unable to create new table";
			exit(0);
		}
	}
	
# compose database population command 			
	open (INPUT,"$datafile") or die "can't open input file $datafile\n $!\n";
	$sql  = " INSERT INTO $tablename VALUES ";
	print STDERR "  Populating table . . .\n";
	
	while (<INPUT>)
	{
		chomp;
	# skip header line
		next if ($_ =~ /taxid/);
		my @fields = split "\t", $_;
		unless (scalar @fields == 2)
		{
			print STDERR "\nSkipping problem line:\n";
			print STDERR "$_\n";
			next;
		}
	# skip taxonomically uninformative entries
		 next unless (defined $fields[1] 
			&& $fields[1] =~ /^\d+$/ 
			&& 	$fields[1] >0
			&& (exists $valid_ids{$fields[1]})						
			);	
		my $full_id = "$fields[0]" || "NULL";
		my $ncbi_tax_id = "$fields[1]"  || "NULL";
		my $species_name = "$valid_ids{$ncbi_tax_id}" || "NULL";			
		
			$sql .= "('$full_id', '$ncbi_tax_id', '$species_name'),\n";
	# send insert statements in pieces to avoid exceeding MySql max_allowed_packet size	
		if ($. % $max_lines_per_packet == 0)
		{
			$sql = substr($sql, 0, -2); # get rid of terminal newline and comma
			$sql .= ";\n";
			$result = $dbh->do($sql);
			$sql  = " INSERT INTO $tablename VALUES ";
		# user feedback for long lists	
			if ($. % 100000 == 0)
			{
				print STDERR "  $. lines processed \n";
			}
		}
		elsif (eof)
		{
			print STDERR "  $. lines processed \n";
			$sql = substr($sql, 0, -2); # get rid of terminal newline and comma
			$sql .= ";\n";
			$result = $dbh->do($sql);
		}
	}	
	
	close INPUT;

# cleanup
	$sth->finish();
	$dbh->disconnect();

###########################
# SUBROUTINES
###########################

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

sub estimate_time
{
	my ($time_mins, $user_message) = @_;
	my $est_time_mins = int($time_mins);
	my $est_time_hours = $time_mins/60;
	$est_time_hours = sprintf "%01X", $est_time_hours;
		
	if ($est_time_hours > 0.99)
	{
		$user_message .= "  NOTE: this step may take $est_time_hours or more hours."; 
	}	
	else
	{
		$user_message .= "  NOTE: this step may take $est_time_mins or more minutes."; 
	}
	return $user_message;
}
__END__