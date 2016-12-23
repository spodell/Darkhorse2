#!/usr/bin/perl
# load_annotation_table.pl
# Sheila Podell
# September 4, 2016

use warnings;
use strict;
use Getopt::Long;
use DBI;

# creates, populates Mysql annotation table based on tab-delimited file in the following format:

# full_id_num   annotation  aa_seq_length 

# will not over-write old table unless -o option specified on command line	
	
# get command line arguments 
	my ($db, $datafile, $sth, $sql, $overwrite, $update, $result, $config_filename);
	my $USAGE;
	my $message = qq(
	Usage: $0  -i source_datafile -c config filename <-o allow overwrites>)."\n\n";
	GetOptions( 
			    "i=s" => \$datafile,
			    "c=s" => \$config_filename,
			    'o'   => \$overwrite,		#allow table overwrites
			    'u'   => \$update,
				);	     
	
	if($USAGE || !$datafile || !$config_filename) 
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

# connect to database
	my $db_path = "DBI:$db_program:$db_name:$db_host:";
	my $dbh = DBI->connect($db_path, $db_user, $db_user_password) or die "Cannot connect: $DBI::errstr";

# check to see if table already exists
	my $tablename = "annotation";
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
	   `full_id` varchar(100) NOT NULL,
	   `annotation` varchar(100)  NULL,
	   `seq_len` int(5) NULL,
		PRIMARY KEY  (`full_id`)
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
	if (defined $update)
	{
		$max_lines_per_packet = 1;
	}

# compose database population command 
	open (INPUT,"$datafile") or die "can't open input file $datafile\n $!\n";
	$sql  = " INSERT INTO $tablename VALUES ";
	print STDERR "  Populating table . . .\n";
	
	while (<INPUT>)
	{
		chomp;
		my @fields = split "\t", $_;
		
	# number of fields should be 3  
		unless (scalar @fields ==3)
		{
			print STDERR "\nSkipping problem line:\n";
			print STDERR "$_\n";
			push @skipped_lines, $_;
			next;
		}
		my $full_id = "$fields[0]";
		my $annotation = "NULL";
		my $seq_length = "$fields[2]";
		 if (defined $fields[1] && length "$fields[1]" > 0)
 		{
 			$annotation = "$fields[1]";
 		}
		else
		{
			print STDERR "no annotation found for $full_id, line $., $_\n";
		}		
		
	# get rid of inconvenient punctuation that chokes MySQL (' and /)
		$annotation =~ s/\'/prime/g;
		$annotation =~ s/\\//g;
				
		$sql .= "('$full_id','$annotation','$seq_length'),\n";
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


__END__