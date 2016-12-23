#!/usr/bin/perl
# load_precalc_lineages.pl
# Sheila Podell
# March 31, 2009

use warnings;
use strict;
use Getopt::Long;
use DBI;

# creates, populates a Mysql table of non-redundant precalculated lineages
# expects input (e.g. from recursive_taxonomy5.pl) as the tab-delimited output
# in the following format:

# species_name ncbi_tax_id	lineage 

# will not over-write old table unless -o option specified on command line	
	
# get command line arguments 
	my ($db, $datafile, $sth, $sql, $overwrite,  $update, $lines_per_packet, $result, $config_filename);
	my $table_name = "precalc_lineages";
	my $USAGE;
	my $message = qq(
	Usage: $0 -d database_name -i source_datafile  -c config filename <-o allow overwrites>)."\n\n";
	GetOptions( "d=s" => \$db,
			    "i=s" => \$datafile,
			    "c=s" => \$config_filename,
			    'o'   => \$overwrite,
			    'u'   => \$update,
				);	     
	
	if($USAGE || !$db || !$datafile || !$config_filename) 
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
	open (INPUT, $config_filename)  or die "can't open DarkHorse config file $config_filename\n$!\n";
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
	my $script_path = $dh_config{script_path} || "not found";
	my $db_name = $dh_config{db_name} || "not found";
	my $db_program = $dh_config{db_program} || "not found";
	my $db_host = $dh_config{db_host} || "not found";
	my $db_user = $dh_config{db_user} || "not found";
	my $db_user_password = $dh_config{db_user_password} || "not found";
	my $max_lines_per_packet = $dh_config{max_lines_per_packet} || "not found";
	my $min_lineage_terms = $dh_config{min_lineage_terms} || "not found";
	
# connect to database
	my $db_path = "DBI:$db_program:$db_name:$db_host:";
	my $dbh = DBI->connect($db_path, $db_user, $db_user_password) or die "Cannot connect: $DBI::errstr";
	
# check to see if table already exists
	my $tablename = "precalc_lineages";
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
			if ($element =~/$table_name/)
			{
				$found_table = 1;
				unless (defined $overwrite || defined $update)
				{
					print STDERR "  Table $table_name already exists.\n";
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
		$sql = qq(CREATE TABLE `precalc_lineages` (
	   `species_name` varchar(100) NOT NULL,
	   `ncbi_tax_id` int(11)  NULL,
	   `lineage` varchar(255)  NULL,
	    `num_terms` int(11) default '0',
	   PRIMARY KEY  (`species_name`),
	   INDEX (`ncbi_tax_id`),
	   INDEX (`num_terms`) 
		) ;\n);
		$result = $dbh->do($sql);
		if ($result)
		{
			print STDERR "New table $table_name created.\n";
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
		$lines_per_packet = 1;
	}

# compose database population command 
	open (INPUT,"$datafile") or die "can't open input file $datafile\n $!\n";
	$sql  = " INSERT INTO $table_name VALUES ";
	print STDERR "  Populating table . . .\n";
	my %uniq_names = ();
	my %uniq_tax_ids = ();
	
	while (<INPUT>)
	{
		chomp;
		next if ($_ =~ /^\s+$/);
		my @fields = split "\t", $_;
		unless (scalar @fields > 0 && scalar @fields < 4)
		{
			print STDERR "\nSkipping problem line:\n";
			print STDERR "$_\n";
			push @skipped_lines, $_;
			next;
		}

		my $species_name = "$fields[0]";
		my $tax_id = "NULL";
		my $lineage = "NULL";
		if (defined $fields[1])
		{
			$tax_id = "$fields[1]";
		}
		if (defined $fields[2])
		{
			$lineage = "$fields[2]";
		}
		my $num_lineage_terms = 0;
		if (defined $lineage)
		{
			my @lineage_terms = split ";", $lineage;
			$num_lineage_terms = scalar @lineage_terms;
			if ($lineage =~ /NULL/)
			{
			
				$num_lineage_terms = 0;
			}
		}
		
	# prevent MySQL failure due to duplicate entries & MySQL case insensitivity
	# some duplicate entries are NCBI taxonomy database errors, 
	# e.g. different tax_id's for same exact species
	# often happens with misspelled species names containing extra underscores, dashes, spaces
		my $uc_species_name = uc($species_name);
		next unless (length $species_name > 2);
		if (exists $uniq_names{$uc_species_name} || exists $uniq_tax_ids{$tax_id})
		{
			#print STDERR "omitting duplicate line\n$_\n";
			next;
		}
		
		$uniq_names{$uc_species_name}++;
		$uniq_tax_ids{$tax_id}++;
		
		$sql .= "('$species_name', '$tax_id', '$lineage', '$num_lineage_terms'),\n";
	# send insert statements in pieces to avoid exceeding MySql max_allowed_packet size	
		if ($. % $lines_per_packet == 0)
		{
			$sql = substr($sql, 0, -2); # get rid of terminal newline and comma
			$sql .= ";\n";
			$result = $dbh->do($sql);
			$sql  = " INSERT INTO $table_name VALUES ";
		# user feedback for long lists	
			if ($. % 50000 == 0)
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
	
# add one final line to take care of cases where no taxonomy ID is available:
	my $species_name = "unknown";
	my $tax_id = 0;
	my $lineage = "NULL";
	my $num_lineage_terms = 0;
	
	$sql = "INSERT INTO $table_name VALUES ";
	$sql .= "('$species_name', '$tax_id', '$lineage', '$num_lineage_terms');\n";
	$result = $dbh->do($sql);

# cleanup
	$sth->finish();
	$dbh->disconnect();

__END__