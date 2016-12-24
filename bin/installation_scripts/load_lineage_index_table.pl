#!/usr/bin/perl
# load_lineage_index_table.pl
# Sheila Podell
# June 3, 2007

# populates a lineage index table based on joining information
# from two pre-existing MySQL tables (db_seq_ids and precalc_lineage).

# This pre-calculated table is needed because the large size of the two 
# existing tables make the join too slow to do as a part of every 
# new DarkHorse run.

# this version of the program takes a table where the join is already done
# by create_missing_lineage_list.pl rather than generating it here

# fields expected in data input:
# 	full_id   gi_num   species_name   lineage   num_terms

# status is "good" or "bad", to indicate the number of lineage terms


use warnings;
use strict;
use Getopt::Long;
use DBI;

# get command line arguments
	my ($db, $datafile, $sth, $sql, $tablename, $found_table, $overwrite, $update,$result, $config_filename);
	my $USAGE;
	my $message = qq(
	Usage: $0 -d database_name -i source_datafile -c config filename <-o allow overwrites>)."\n\n";
	GetOptions( "d=s" => \$db,
			    "i=s" => \$datafile,	
			    "c=s" => \$config_filename,
			    'o'   => \$overwrite,		
			    'u'   => \$update,
				);	     
	
	if($USAGE || !$db) 
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
	my $path = "DBI:mysql:$db:vitro.ucsd.edu:";
	my $db_path = "DBI:$db_program:$db_name:$db_host:";
	my $dbh = DBI->connect($db_path, $db_user, $db_user_password) or die "Cannot connect: $DBI::errstr";
	my @row = ();
	my $element;
	
# check to see if lineaage index table already exists
	$tablename = "lineage_index";
	$found_table = 0;
	$sql = "show tables;";
	$sth = $dbh->prepare($sql)
    	or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
	$sth->execute()
    	or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");
	@row = ();
	$element = "";
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
					#print STDERR "  Specify -o option if you wish to overwrite.\n";
					#print STDERR "  Specify -u option to update without overwriting.\n";
					exit(0);
				}
				last;
			}
		}
	}
	
# open the tab file, use to create new table
	unless (-s $datafile)
	{
		print STDERR "  Couldn't find source datafile $datafile \n";
		exit(0);
	}	

# if over-write option specified, drop table from database.

if (defined $overwrite || $found_table == 0)
	{
		unless ($found_table == 0)
		{
			print STDERR "Over-writing table $tablename.\n";
		}
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
       `full_id` varchar(100) NOT NULL,
       `ncbi_tax_id` int(11)  NULL,
       `species_name` varchar(100) NULL,
       `lineage` varchar(255) NULL,
       `num_terms` int(3) NOT NULL,
        PRIMARY KEY  (`full_id`),
        INDEX (`num_terms`)
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

# figure out how many lines need to be entered, allowing insert statements to 
# be sent in pieces (avoid exceeding MySql max_allowed_packet size)
	my $datafile_wc = `wc -l $datafile`;
	$datafile_wc  =~ /(\d+).+/;
	my $total_num = $1;
	print STDERR "total num = $total_num\n";
	my @skipped_lines = ();
	if (defined $overwrite || $found_table == 0)
	{
		#$max_lines_per_packet = 2000;
	}
	elsif (defined $update)
	{
		$max_lines_per_packet = 1;
	}

# compose database population command, send insert statement packets
	open (INPUT,"$datafile") or die "can't open input file $datafile\n $!\n";
	$sql  = " INSERT INTO $tablename VALUES ";
	print STDERR "  Populating table . . .\n";
	
	while (<INPUT>)
	{
		chomp;
		my @fields = split "\t", $_;
		unless (scalar @fields == 5)
		{
			print STDERR "\nSkipping problem line:\n";
			print STDERR "$_\n";
			exit(0);
			push @skipped_lines, $_;
			next;
		}
		my $full_id = "$fields[0]";
		my $ncbi_tax_id = "$fields[1]";
		my $species_name = "$fields[2]";
		my $lineage = "$fields[3]";
		my $num_terms = "$fields[4]";
		
	# send insert statements in pieces to avoid exceeding MySql max_allowed_packet size	
		$sql .= "('$full_id','$ncbi_tax_id','$species_name','$lineage', '$num_terms'),\n";
		if ($. % $max_lines_per_packet == 0)
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


__END__
