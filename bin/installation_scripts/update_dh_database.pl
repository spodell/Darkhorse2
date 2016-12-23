#!/usr/bin/perl
# update_dh_database.pl
# Sheila Podell 
# December 19, 2016

# updates DarkHorse reference database with new sequences  
# takes as input 
	# DarkHorse config file, specifying pre-existing database
	# tab-delimited file containing the following columns:
		# protein_id_num
		# Species name
		# colon-delimited lineage string from superkingdom to genus level, e.g.
			# Bacteria;Cyanobacteria;Gloeobacteria;Gloeobacterales;Gloeobacteraceae;Gloeobacter
	# fasta format amino acid sequence file with sequence ids matching tab-delimited file

# output:
	# modified mysql database
	# new fasta file with taxonomically informative ids

use warnings; 
use strict;
use Getopt::Long;
use Cwd;
use File::Basename;
use File::Spec;
use DBI;

# global input:
    my ($db, $sth, $sql, $overwrite, $update, $result, $msg);
    my ($table_name, $firstline, $datafile, $fastafile);
    my $config_template = "templates/config_template";
    my $config_filename = "$config_template";
    my $user_message = "";
    my $cmd = "";
    my $debug = 0;
        
    my $USAGE;
	my $message = qq(
  Usage: $0 -i input datafile filename -f fasta filename -c configuration file\n\n); 
    
	GetOptions( 
			    "i=s" => \$datafile,
			    "f=s" => \$fastafile,
			    "c=s" => \$config_filename,   # pre-existing or manually edited file
			    'o'   => \$overwrite,		
			    'u'   => \$update,
			    "d=i" => \$debug,
				);
	 &check_args();
	 
# get SQL database parameters from config file
    my %config_params= &read_config_file($config_filename);
    my $program_directory = $config_params{program_directory};
	my $db_name = $config_params{db_name} || "not found";
	my $db_program = $config_params{db_program} || "not found";
	my $db_host = $config_params{db_host} || "not found";
	my $db_user = $config_params{db_user} || "not found";
	my $db_user_password = $config_params{db_user_password} || "not found";
	my $max_lines_per_packet = $config_params{max_lines_per_packet} || "not found";
			
# connect to database
	my $db_path = "DBI:$db_program:$db_name:$db_host:";
	my $dbh = DBI->connect($db_path, $db_user, $db_user_password) or die "Cannot connect: $DBI::errstr"; 

# check to make sure that required database tables are already in place
	my @table_list = (
	"db_seq_ids",
	"lineage_index",
	"precalc_lineages",);
	foreach my $tbl (@table_list)
	{
		my $found_table = &check_table($tbl);
		unless ($found_table ==1)
		{
			print STDERR "ERROR: could not find required table $tbl in database $db_name\n";
			exit(0);
		}		
	}
		
# create logfile to keep track of progress
	# my $logfile = "darkhorse_update_"."$$".".log";
# 	open (LOGFILE, ">$logfile") or die "couldn't open logfile $logfile for writing, $!\n";	
# 	print LOGFILE qq(###################################################
#   Darkhorse version 2.0 (beta) update logfile
# ###################################################
# 
# Installation settings:
# 	config filename = $config_filename
# 	MySQL database name = $config_params{db_name}
# 	max_lines_per_packet = $config_params{max_lines_per_packet} 
# 
# );


# Get lineages and corresponding tax ids already in database (precalc_lineages table)
	my %lineages = &get_precalc_lineages;	# key = lineage, value = tax_id

# get protein information from fasta file, put into hash table
	my $fastafile_size = &check_filesize($fastafile); # MB
	my $filesize_txt = sprintf("%.1f", $fastafile_size);
	my $trouble_flag = $fastafile_size/$max_lines_per_packet;	
	if ($trouble_flag > 20)
	{
		$msg = "    WARNING: Fasta filesize is $filesize_txt MB. Might need to split into smaller pieces to avoid memory errors.";
		&send_message($msg);
	}
	 
	my %fasta_info = &get_prot_info($fastafile);	# key = id number  # value = sequence_object
	
# create fasta output file
	my $fasta_out = "informative_update_sequences.faa";
	open (FASTA, ">$fasta_out") or die "unable to open fasta output file $fasta_out for writing,$!\n";
	
# parse input data file	
	my %species_lineage_index = (); # key = species name, value = lineage
	my %name_count = ();
	my %id_nums = ();	#key = id num # value # times seen (eliminate duplicates)	
	my $fake_tax_id = 0;
	my ($current_tax_id, $num_lineage_terms);
	my $num_bad = 0;
	my $num_db_seq_ids_entries = 0;
	my $num_lineage_index_entries = 0;
	my $num_new_lineages = 0;
				
	open (INFILE, "$datafile") or die "couldn't open $datafile, $!\n";	 
	while (<INFILE>)
	{		
		if ($num_bad > 20)
		{
			$msg = "    PROGRAM EXIT: Too many errors.\n";
			&send_message($msg);
			exit(0);
		}
		chomp;
		next if ($_ =~ /^\s+$/);
		my @dataline = split '\t', $_;
		
		unless (scalar @dataline == 3)
		{
			$msg = "    WARNING: skipping datafile line $.: incorrect number of columns.";
			$num_bad++;
			&send_message("$msg");
			next;
		}
		
		my $current_id_num = $dataline[0];
		my $current_name = $dataline[1];
		my $current_lineage = $dataline[2];	
	
	# be sure the id number entry is not a duplicate
		$id_nums{$current_id_num}++;
		if ($id_nums{$current_id_num} > 1)
		{
			$msg = "    WARNING: skipping datafile line $. (ID_NUM = $current_id_num). Duplicate id number.";
			&send_message("$msg");
			next;
		}
	
	# get rid of un-necessary extra spaces
		$current_lineage =~ s/; /;/g;
		$current_lineage =~ s/ ;/;/g;
			
	# check to be sure there was a fasta protein file
		unless (exists $fasta_info{$current_id_num})
		{
			$msg = "    WARNING: skipping datafile line $. (ID_NUM = $current_id_num). No fasta protein file found.";
			$num_bad++;
			&send_message($msg);
			next;
		}	
	
	# check validity of current lineage, get rid of uninformative terms
	# skip with error message if format is bad
		$current_lineage = &check_lineage_format($current_lineage);
		unless (defined $current_lineage && length $current_lineage > 2)
		{
			$msg = "    WARNING: skipping datafile line $.: Invalid lineage format.\n\t$current_lineage";
			&send_message($msg);
			$num_bad++;
			next;
		}
		
	# get number of terms in revised lineage (uninformative terms removed)
			my @current_terms = split ';', $current_lineage;
			$num_lineage_terms = scalar @current_terms;
	
	# check for uniqueness of name/lineage pair 
		if (exists $species_lineage_index{$current_name})
		{
			unless ("$current_lineage" eq "$species_lineage_index{$current_name}")
			{
			 	$msg = "    FATAL ERROR: inconsistent lineage, datafile line $. . 
\tID_NUM = $current_id_num 
\tNAME = $current_name 
\tCURRENT_LINEAGE = $current_lineage
\tPREVIOUS_LINEAGE = $species_lineage_index{$current_name}";
			 	&send_message($msg);
			 	exit(0);
			}
		}
		else
		{
			$species_lineage_index{$current_name} = $current_lineage;
		}
		
	# first time name/lineage pair is seen, check to see if they're already in the database
		$name_count{$current_name}++;
		if ($name_count{$current_name} < 2)  
		{		
			if (exists $lineages{$current_lineage})
			{	
				# use same tax_id as lineage already in database
				# even though may not be isn't strictly correct, it won't affect DarkHorse program function
				# and will take care of cases where no ncbi_tax_id has yet been assigned
				$current_tax_id = $lineages{$current_lineage};
				
				# if necessary species name isn't already in precalc_lineages table, add it				
				$result = &check_precalc_lineage_species_name($current_name);
				unless ($result > 0)
				{
					$msg = "    Adding new species name to precalc lineages table: $current_name ";
					&send_message($msg);
					$sql  = " INSERT INTO precalc_lineages VALUES ";
					$sql .= "('$current_name', '$current_tax_id', '$current_lineage', '$num_lineage_terms');";					
					$result = &sql_execute($sql);
					if ($result >0)
					{
						$num_new_lineages++;
					}
					else
					{
						$num_bad++;
						next;
					}
					$result = 0;
				}
			}
			else
			{
			# if lineage isn't already there, need to invent new tax id for it (negative number)
			# make these successively more negative for additional new lineages
				$current_tax_id = $fake_tax_id -1;
				$lineages{$current_lineage} = $current_tax_id;
				$fake_tax_id = $current_tax_id;

			# here, need to update precalc_lineages table with both name and lineage				
				$msg = "    WARNING: Lineage not found in current database. Creating placeholder for tax_id = $fake_tax_id";
				$msg .= "    $current_lineage";
				&send_message($msg);
						
				my $tax_id = "$lineages{$current_lineage}";
				    	
			 
				$sql  = " INSERT INTO precalc_lineages VALUES ";
				$sql .= "('$current_name', '$tax_id', '$current_lineage', '$num_lineage_terms');";
				$result = &sql_execute($sql);
				if ($result >0)
				{
					$num_new_lineages++;
				}
				else
				{
					$num_bad++;
					next;
				}
				$result = 0;
				$sql = "";						
			}		
		}
	# update db_seq_ids	table
		my $annotation = "$fasta_info{$current_id_num}->{annotation}"; #my %fasta_info = &get_prot_info($fastafile);	# key = id number  # value = sequence_object
		my $seq_length = "$fasta_info{$current_id_num}->{seq_length}";
		
		$sql  = " INSERT INTO db_seq_ids VALUES ";
		$sql .= "('$current_id_num','$current_tax_id','$current_name','$annotation','$seq_length');\n";
		#print STDERR "$sql\n";
		$result = &sql_execute($sql);
		if ($result >0)
		{
			$num_db_seq_ids_entries++;
		}
		else
		{
			$num_bad++;
			next;
		}
		$result = 0;
		$sql = "";
		
	# update lineage_index table
		$sql  = " INSERT INTO lineage_index VALUES ";
		$sql .= "('$current_id_num','$current_tax_id','$current_name','$current_lineage', '$num_lineage_terms');\n";	
		$result = &sql_execute($sql);
		if ($result >0)
		{
			$num_lineage_index_entries++;
		}
		else
		{
			$num_bad++;
			next;
		}
		$result = 0;
		$sql = "";
		
	# write successfully updated sequences to fasta output file  
		my $outseq = &pretty_seq("$fasta_info{$current_id_num}->{sequence}");
		my $sequence_txt = ">"."$current_id_num"." $annotation\n"."$outseq";
		print FASTA "$sequence_txt\n";										
	}
		
# user feedback on number of records added to each table
	$msg = "    Finished $db_name updates.\n";
	$msg .= "      $num_new_lineages new entries added to precalc_lineages table.\n";
	$msg .= "      $num_db_seq_ids_entries new entries added to db_seq_ids and lineage_index tables.";		
	&send_message($msg);

# clean up: close/remove temp files, close database connections
	close (FASTA);	
	 
#########################################################
sub check_args
{
 	if($USAGE || !$datafile || !$config_filename || !$fastafile) 
	{
		print STDERR $message;
		exit(0);
	} 
	unless (-s "$datafile")
	{
		print STDERR "  Couldn't find input data file $datafile\n";
		exit(0);
	}

	unless (-s "$config_filename")
	{
		print STDERR "  Couldn't find config file $config_filename \n";
		exit(0);		
	}
	
	unless (-s "$fastafile")
	{
		print STDERR "  Couldn't find fastafile $fastafile \n";
		exit(0);		
	}	 
}

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

sub check_db_connection
{
    my ($hashref) = @_;
	my %dh_config = %$hashref;
	my $found_db = 0;
	
	my $db_name = $dh_config{db_name} || "not found";
	my $db_program = $dh_config{db_program} || "not found";
	my $db_host = $dh_config{db_host} || "not found";
	my $db_user = $dh_config{db_user} || "not found";
	my $db_user_password = $dh_config{db_user_password} || "not found";
    
    my $db_path = "DBI:$db_program:$db_name:$db_host:";
    
	my $dbh = DBI->connect($db_path, $db_user, $db_user_password);
	unless ($dbh)
	{
	    return $found_db;	    
	}

	my $sql = "show databases";
	my $sth = $dbh->prepare($sql)
    	or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
	$sth->execute()
    	or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");
	my @row = ();
	my $element = ();
	while (@row = $sth->fetchrow_array)
	{
		foreach $element (@row)
		{
			if ($element =~/$db_name/)
			{
				$found_db = 1;
				last;
			}
		}
	}	
	$dbh->disconnect();
	return $found_db;
}

sub check_table
{
	my ($tablename) = @_;
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
				last;
			}
		}
	}
	return $found_table;
}

#check_precalc_lineage_species_name($current_name);
sub check_precalc_lineage_species_name
{
	my ($name) = @_;
	my $found = 0;
	$sql = "select species_name from precalc_lineages where species_name = '$name';";
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
			if ($element =~/$name/)
			{
				$found = 1;
				last;
			}
		}
	}
	return $found;
}

sub sql_execute
{
	my ($sql_statement) = @_;
	my $sql_result = $dbh->do($sql_statement);
	my $message = "";
	my $success = 0;
	
	if (defined $sql_result && $sql_result > 0)
	{
		$message = "     SUCCESSFUL SQL: $sql_statement";
		$success++;
	}
	else
	{
		$message = "     SQL FAILURE: $sql_statement";
		&send_message($message);
	}
	return($success);

}

sub send_message
{
	my ($msg, $error) = @_;
	
#	print LOGFILE "$msg\n";
	print STDERR "$msg\n";	
}

sub new_seq {
  my ($className, $param) = @_;
  my $self = {};
  bless $self, $className;
  my @properties = @$param;
 
  $properties[0] =~ s/\>//;
  $self->{id} = shift @properties;   
  $self->{annotation} = join " ", @properties;
  
  # get rid of inconvenient punctuation 
  	$self->{annotation} =~ s/\'//g;
	$self->{annotation} =~ s/\\//g;
	$self->{annotation} =~ s/\;/,/g;
	$self->{annotation} =~ s/\"//g;
	$self->{annotation} =~ s/RecName: Full=//;
  
  $self->{sequence} ="";
  $self->{seq_length} =""; 
  
  return($self)
}
	
sub check_lineage_format
{
	my ($lineage) = @_;
	my $validity = 0;
	$lineage =~ s/; /;/g;
	my @terms = split ";", $lineage;
	
	my @informative_terms = ();
	foreach my $next_term (@terms)
	{
		$next_term =~ s/unclassified//;
		$next_term =~ s/\(miscellaneous\)//;
		next if  $next_term =~ /candidate\s/;
		next if  $next_term =~ /\(class\)/;
		next if  $next_term =~ /group/;
		next if  $next_term =~ /subdivision/;
		next if  $next_term =~ /subgen/;
		next if  $next_term =~ /subg\./;
		next if  $next_term =~ /Family/;
		
		 if (defined $next_term 
		&& length $next_term > 2
		&& $next_term =~ /^\w+/)
		{
			push @informative_terms, $next_term;
		}	
	}
	
	my $num_good = scalar @informative_terms;
	unless ($num_good > 0)
	{
		return(0);
	}
	my $revised_lineage = join ";", @informative_terms;;
	return $revised_lineage;	
}

sub get_prot_info
{
	my ($infilename) = @_;
	my %seq_objects = ();
	
	open (INFILE, "<$infilename") or die "couldn't open $infilename, $!\n";
	my $current = "";
	
	while (my $line = <INFILE>)
	{
		chomp $line;
		if (eof) 
		{					
			unless ($line =~ /^\s+$/)
			{
				$current->{sequence} .= uc"$line";
			}			
			$current->{seq_length} = length "$current->{sequence}";
			my $id = $current->{id};
						
			my $protein_check = &check_sequence_type("$current->{sequence}");
			unless ($protein_check eq "protein")
			{
				$msg = "    WARNING: skipping fasta file $id. Fasta protein sequence format error.";
				&send_message($msg);
				next;
			}
			
			$seq_objects{$id} = $current;
				last;
		}
		
		next if $line =~ /^\s+$/;		
		if ($line =~ />(.+)/ || eof) # this prints last fasta at end of file without blank line
		{								
		# process previous object unless it hasn't been created yet (1st line of file)
			unless ($. < 2)
			{				
				my $id = $current->{id};

		# check to be sure it's a protein
			my $protein_check = &check_sequence_type("$current->{sequence}");
			unless ($protein_check eq "protein")
			{
				$msg = "    WARNING skipping id number $id. Fasta protein sequence error.";
				&send_message($msg);
				next;
			}
			$current->{seq_length} = length "$current->{sequence}";
			$seq_objects{$id} = $current;				
			}
			
			my @tmp = split " ", $line;
			$current = new_seq("record",\@tmp);
		}													
		else
		{			
			$current->{sequence} .= uc"$line";
		}
	}
	
	
	close INFILE;
	return %seq_objects;		
}
	
sub check_filesize
{
	my ($name) = @_;
	my $filesize = -s $name;	
	my $filesize_KB = $filesize/1000;
	my $filesize_MB = $filesize_KB/1000;
	my $filesize_GB = $filesize_MB/1000;
	my $filesize_txt = "$filesize_KB";
		$filesize_txt .= " KB";
	if ($filesize_MB >1)
	{
		$filesize_txt = sprintf("%.1f", $filesize_MB);
		$filesize_txt .= " MB";
	}
	if ($filesize_GB >1)
	{
		$filesize_txt = sprintf("%.1f", $filesize_GB);
		$filesize_txt .= " GB";
	}
	
	return $filesize_MB;
}

sub check_sequence_type
{
	my ($seqstring) = @_;
	$seqstring =~ s/\s//g;
	return unless (defined $seqstring && length $seqstring > 2);
	
	my $seqtype = "";
		
	my $seqlength = length $seqstring;	
	my $num_ATCGN = (uc$seqstring =~ tr/ATUCGN//);
	my $num_aminos = (uc$seqstring =~ tr/EFILPQ//) || 0;
	
	if  ($num_ATCGN eq  $seqlength || $num_aminos < 1)
	{
		$seqtype = "nucleic";
		return $seqtype;
	}			
	elsif ($num_aminos > 1)
	{
		$seqtype = "protein";
		return $seqtype;
	}	
	else
	{
		$seqtype = "ambiguous";
		return $seqtype;
	}
	
	return $seqtype;
}

sub get_precalc_lineages
{
	my %lineage_hash = ();
	my $tablename = "precalc_lineages";
	$sql = "select lineage, ncbi_tax_id from precalc_lineages;";
	$sth = $dbh->prepare($sql)
    	or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
	$sth->execute()
    	or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");
	my @row = ();
	my $element = ();
	while (@row = $sth->fetchrow_array)
	{
		$lineage_hash{$row[0]} = $row[1];
	}
	return %lineage_hash;
}

sub pretty_seq
{
	my ($seq) = @_;
	my $pretty_seq = "";
	my $count = 0;
	my $line_length = 70;
	$seq =~ s/ //g;
	while ($count < length $seq)
	{
		$pretty_seq .= substr ($seq, $count, $line_length);
		$pretty_seq .= "\n";
		$count += $line_length;
	}
	chomp $pretty_seq;
	return $pretty_seq;
}



