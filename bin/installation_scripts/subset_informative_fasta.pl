#!/usr/bin/perl
# subset_informative_fasta.pl
# Sheila Podell
# December 23, 2016

# using DarkHorse MySQL tables and matching taxonomically informative fasta sequences
# selects those sequences whose lineages match a list of input keywords

# requires DarkHorse accessory script
#	bin/installation_scripts/getseq_multiple_nr.pl

# outputs the selected sequences in fasta format to output file

use warnings;
use strict;
use Getopt::Long;
use DBI;
use Cwd;
use File::Basename;
use File::Spec;

# get command line arguments
	my ($db, $sth, $sql, $tablename, $found_table, $overwrite, $update, $result);
	my ($fastafile, $config_filename, $keywords_file);
	my $USAGE;
	my $message = qq(
	Usage: $0  -i keywords_filename -c config_filename -f source_fasta_filename)."\n\n";
	GetOptions( 
			    "i=s" => \$keywords_file,	
			    "c=s" => \$config_filename,
			    "f=s" => \$fastafile,
				);	     	
	&check_args();
	
# get SQL database parameters from config file
    my %config_params= &read_config_file($config_filename);
    my $program_directory = $config_params{program_directory} || "not found";
	my $db_name = $config_params{db_name} || "not found";
	my $db_program = $config_params{db_program} || "not found";
	my $db_host = $config_params{db_host} || "not found";
	my $db_user = $config_params{db_user} || "not found";
	my $db_user_password = $config_params{db_user_password} || "not found";
	my $max_lines_per_packet = $config_params{max_lines_per_packet} || "not found";
	my $min_lineage_terms = $config_params{min_lineage_terms} || "not found";
	
# convert relative to absolute paths
	my $abs_prog_directory = File::Spec->rel2abs($program_directory);
	my $abs_config_file = File::Spec->rel2abs($config_filename);
	my $abs_fastafile =  File::Spec->rel2abs($fastafile);
	my $abs_keywords_file = File::Spec->rel2abs($keywords_file);

# check for required accessory program 
	my $getseq_program = ("$abs_prog_directory/bin/installation_scripts/getseq_multiple_nr.pl");
	unless (-s "$getseq_program" && -x "$getseq_program")
	{
		print STDERR "Couldn't find required program $getseq_program\n";
		exit (0);
	}         
	
# connect to database
	my $db_path = "DBI:$db_program:$db_name:$db_host:";
	my $dbh = DBI->connect($db_path, $db_user, $db_user_password) or die "Cannot connect: $DBI::errstr";
	my @row = ();
	my $element;

# check to make sure lineage index table already exists
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
				last;
			}
		}
	}
	
	unless ($found_table > 0)
	{
		print STDERR "  FATAL ERROR: unable to find table lineage_index in database $db_name\n";
		exit (0);
	}
	
# get list of keywords from input file, put into tmp output
	my @sql_search_terms = ();
	open (INPUT, "$abs_keywords_file") or die "unable to open keywords file $abs_keywords_file.\n$!\n";
	while (<INPUT>)
	{
		$_ =~ s/\r\n/\n/s; # convert windows to unix line endings, if necessary
		chomp;
		next if ($_ =~ /^\s+$/ || length $_ < 2);
		my $wildcard_kw = qq("%$_%");
		push @sql_search_terms, "$wildcard_kw";
	}
	close INPUT;
	
# get list of id numbers by sql query of lineage_index table from database into tmp file
	my $tmp_id_list_filename = "tmp_id_list";
	my $num_selected = &select_sequence_ids($tmp_id_list_filename, \@sql_search_terms);
	print STDERR "  Found $num_selected id numbers matching criteria\n";

# make sure id numbers in tmp id list file are unique
	my $tmp2 = "$tmp_id_list_filename"."_uniq";
	`sort $tmp_id_list_filename |uniq > $tmp2`;
	`mv $tmp2 $tmp_id_list_filename`;
	my $uniq_hits = `wc -l $tmp_id_list_filename`;
	chomp $uniq_hits;
	
	print STDERR "  Found $uniq_hits unique/$num_selected id numbers matching criteria\n";	

# select sequences using accessory script, send to file 
	my $select_fasta_filename = "tax_select_dh.fasta";
	my $cmd = "$getseq_program $tmp_id_list_filename $abs_fastafile > $select_fasta_filename\n";
	print STDERR "  selecting $num_selected fasta sequences\n";
	`$cmd`;
	unlink $tmp_id_list_filename;
	
#########################################################
# SUBROUTINES
#########################################################
sub check_args
{
 
	if ($USAGE || !$config_filename || !$keywords_file || !$fastafile) 
	{
		print STDERR $message;
		exit(0);
	}
	unless (-s $config_filename)
	{
		print STDERR "  Couldn't find config file $config_filename \n";
		exit(0);
	}
	unless (-s $keywords_file)
	{
		print STDERR "  Couldn't find list of input taxonomy terms $keywords_file \n";
		exit(0);
	}
	unless (-s $fastafile)
	{
		print STDERR "  Couldn't find source fasta file $fastafile \n";
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

sub select_sequence_ids
{
	my ($outfile, $list_arrayref) = @_;
	my @kw_list = @$list_arrayref;
	my $num_hits = 0;
	
	# results too big to store hash (too much RAM) - send direct to file instead
	
	open (OUTPUT, ">$outfile") or die "can't open temp file $outfile for writing\n $!\n";

	# set up db connection;
		my $db_name = $config_params{db_name};
		my $db_program = $config_params{db_program};
		my $db_host = $config_params{db_host};
		my $db_user = $config_params{db_user};
		my $db_user_password = $config_params{db_user_password};    
		my $db_path = "DBI:$db_program:$db_name:$db_host:";
		my $dbh = DBI->connect($db_path, $db_user, $db_user_password);

	# run separate query for each keyword, append to output file	
	foreach my $kw (@kw_list)
	{
		$sql  = "select full_id from lineage_index where lineage like $kw;";

		print STDERR "  executing query: $sql\n";
		$sth = $dbh->prepare($sql)
		or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
		$sth->execute()
		or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");
	
		my @row = ();
		my $element = ();
		
		while (@row = $sth->fetchrow_array)
		{									
			my $current_id = $row[0];
			next unless (defined $current_id && length $current_id >0);
			print OUTPUT "$current_id\n";
			$num_hits++;
		}
	}

	$dbh->disconnect();
	close OUTPUT;
	return $num_hits;
}




	
	
