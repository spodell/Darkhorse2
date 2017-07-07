#!/usr/bin/env perl
# load_taxonomy_names.pl
# Sheila Podell
# June 19, 2007

use warnings;
use strict;
use Getopt::Long;
use DBI;

# creates, populates a relational database table of taxonomy names,
# based on NCBI Taxonomy database file "names.dmp", downloaded from NCBI:

	# ncftpget ftp://ftp.ncbi.nih.gov/../pub/taxonomy/taxdump.tar.gz 

# will not over-write old table unless -o option specified on command line	
	
# get command line arguments, defaults 
	my ($db, $datafile, $sth, $sql, $overwrite, $update, $lines_per_packet, $result, $config_filename);
	my $USAGE;
	my $message = qq(
	Usage: $0 -d database_name -i source_datafile -c config filename <-o allow overwrites> )."\n\n";
	GetOptions( "d=s" => \$db,
			    "i=s" => \$datafile,
			    "c=s" => \$config_filename,
			    'o'   => \$overwrite,
			    'u'   => \$update,
				);	     
	
	if($USAGE || !$db || !$datafile) 
	{
		print STDERR $message;
		exit(0);
	} 
	unless (-s $datafile)
	{
		print STDERR "  Couldn't find source datafile $datafile \n";
		exit(0);
	}
	unless (-s $config_filename)
	{
		print STDERR "  Couldn't find config file $config_filename \n";
		exit(0);
	}

# get SQL database parameters from config file
    my %config_params= &read_config_file($config_filename);
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
	
# check to see if table already exists
	my $tablename = "names";
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
		print STDERR "\nOver-writing table $tablename.\n";
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
		$sql = qq(CREATE TABLE `names` (
		`id` int(11) NOT NULL auto_increment,
		`tax_id` int(11) NOT NULL default '0',
		`name_txt` varchar(255) default NULL,
		`unique_name` varchar(255) default NULL,
		`name_class` varchar(255) default NULL,
		PRIMARY KEY  (`id`),
		KEY `tx` (`tax_id`),
		KEY `name` (`name_txt`)
);\n);
			
		
		$result = $dbh->do($sql);
		if ($result)
		{
			print STDERR "New table $tablename created.\n";
		} 
		else
		{
			print STDERR "Unable to create new table";
			exit(0);
		}
    }
    
# figure out how many lines need to be entered, allowing insert statements to 
# be sent in pieces (avoid exceeding MySql max_allowed_packet size)
	my $datafile_wc = `wc -l $datafile`;
	$datafile_wc  =~ /(\d+).+/;
	my $total_num = $1;
	
	print STDERR "total num = $total_num\n";
	my @skipped_lines = ();
	if (defined $overwrite || $found_table == 0)
	{
		$lines_per_packet = 2000;
	}
	elsif (defined $update)
	{
		$lines_per_packet = 1;
	}

# compose database population command 
	open (INPUT,"$datafile") or die "can't open input file $datafile\n $!\n";
	$sql  = " INSERT INTO $tablename VALUES ";
	print STDERR "  Populating table . . .\n";
	my %uniq_lines = ();
	
	while (<INPUT>)
	{
		chomp;
		next if ($_ =~ /^\s+$/);
	# get rid of extra single quotes and tab characters
		$_ =~ s/'//g;
		$_ =~ s/\t//g;
			
		my @fields = split /\|/, $_;
		
		unless (scalar @fields > 2 )
		{
			print STDERR "\nSkipping problem line:\n";
			print STDERR "$_\n";
			push @skipped_lines, $_;
			next;
		}

		my $tax_id = "$fields[0]";
		my $name_txt = "$fields[1]";
		my $unique_name = "$fields[2]";
		my $name_class = "$fields[3]";
	
	# first field blank to accomodate auto-increment	
		$sql .= "('0', '$tax_id', '$name_txt', '$unique_name', '$name_class'),\n";
	
	# send insert statements in pieces to avoid exceeding MySql max_allowed_packet size	
		if ($. % $lines_per_packet == 0)
		{
			$sql = substr($sql, 0, -2); # get rid of terminal newline and comma
			$sql .= ";\n";
			$result = $dbh->do($sql);
			$sql  = " INSERT INTO $tablename VALUES ";
		# user feedback for long lists	
			if ($. % 50000 == 0)
			{
				print STDERR "  $. lines processed \n";
			}
		}
		elsif ($. == $total_num)
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
__END__
