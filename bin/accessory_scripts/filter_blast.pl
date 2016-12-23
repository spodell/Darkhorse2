#!/usr/bin/perl
# filter_blast.pl 
# Sheila Podell 
# July 19,2010

# Takes as input 3 files:
#	fasta query file used for blast
#	darkhorse.cfg file, containing darkhorse database access information
# 	(to obtain pre-calculated lengths for match sequences) 
#	blast output (tab-delimited) m8 format

# Finds amino acid lengths for each query-match pair
# Removes pairs with alignment length less than specified percentage
# of either query or match 

# requires MySQL database containing information on match sequences
# (table gb_seq_ids should be loaded as part of DarkHorse program installation)

# outputs revised file of blast hits containing only those rows
# where alignment length is > specified percentage of both query and match 

# Also outputs two tab-delimited flat files obtained during database queries
# saved for future DarkHorse use, to avoid later repetition  
	# query_info: full_id   protein   size   annotation (if avail.)
	# match_info: full_id  tax_id  species_nam  size 

use warnings; 
use strict;
use Getopt::Long;
use DBI;
use Benchmark;
use File::Basename;

# set default variables
	my $t0 = new Benchmark;
	my $query_fasta = "";
	my $blast_tab = "";
	my ($db_name, $db_program, $db_host, $db_user, 
		$db_user_password, $sth, $max_lines_per_packet, $min_lineage_terms,
		$pct_cutoff, $config_filename) ;

# Get command line options
	my $USAGE;
	my $message = qq(
  Usage: $0 -q query_fasta_file -t blast_tab_input -c config_file
	  
  Optional parameters (override config_file):
	-f  <filter threshold for minimum alignment length coverage> 	    
)."\n";

# get command line arguments
	GetOptions( "q=s" => \$query_fasta,		# INPUT name					
				"t=s" => \$blast_tab,		# tab-delimited (m8) blast filename
				"c=s" => \$config_filename,		# tab-delimited (m8) blast filename
				"f=f" => \$pct_cutoff,		# minimum alignment coverage for keepers					
        );
   
 # get global parameters from config file
	unless (defined $config_filename)
		{
			print $message;
		    exit(0);
		}
	 unless (-s $config_filename)
     {
        print STDERR "Couldn't read config file $config_filename\n";
        exit(0);
     }
	
	 my %dh_config = &read_config_file("$config_filename");
	 
	$db_name = $dh_config{db_name} || "not found";
	$db_program = $dh_config{db_program} || "not found";
	$db_host = $dh_config{db_host} || "not found";
	$db_user = $dh_config{db_user} || "not found";
	$db_user_password = $dh_config{db_user_password} || "not found";
	$max_lines_per_packet = $dh_config{max_lines_per_packet} || "not found";
	unless (defined $pct_cutoff)
	{
		$pct_cutoff =$dh_config{min_align_coverage};
	}
	
	& check_args(); 
	
# connect to database
	my $db_path = "DBI:$db_program:$db_name:$db_host:";
	my $dbh = DBI->connect($db_path, $db_user, $db_user_password) or die "Cannot connect: $DBI::errstr";

# get and store match info from darkhorse database
	my $t_start_sql_query = new Benchmark;
	print STDERR "getting match sequence lengths\n";
	my %match_lengths = ();	#key = full_id, value = $ncbi_tax_id\t$species_name\t$seq_length\t$full_id	
	my $getlist = `cut -f2 $blast_tab | cut -f2 -d '|' |sort| uniq`;	
	my @getlist = split "\n", $getlist;	
	
	my $table_name = "db_seq_ids";
	my $firstline = "SELECT full_id, ncbi_tax_id, species_name, seq_length FROM $table_name where full_id in (";		
	my @species_table=&sql_select(\@getlist,$table_name, $firstline);
	
	my $t_finish_sql_query = new Benchmark;
	
	my $t_start_hash_convert = new Benchmark;
	
	foreach (@species_table)
	{		
	#skip header lines
		next if ($_ =~ /^\s/);
		next if $_ =~ /query_id/;
		next if $_ =~ /query/i;
			
		my ($full_id, $ncbi_tax_id, $species_name, $seq_length, $gi_num) = split "\t", $_;
		$match_lengths{$full_id} =  "$ncbi_tax_id\t$species_name\t$seq_length\t$full_id";		
	}
	my $t_finish_hash_convert = new Benchmark;
	my $num_matchlengths = scalar (keys %match_lengths);
	 
# parse query file sequences
	my $t_start_get_query_fasta = new Benchmark;
	print STDERR "getting query fasta sequences\n";
	open (INPUT, "<$query_fasta") or die "couldn't open $query_fasta, $!\n";
	
	my $flag = 0;
	my %seqs = (); #list of sequence objects key = id value = object reference
	my $current = "";
	while (my $line = <INPUT>)
	{
		next if $line =~ /^\s+$/;
		chomp $line;
		
		if ($line =~ />(.+)/) 
		{								
			my @tmp = split " ", $line;
			$current = new_seq("record",\@tmp);
			$seqs{$current->{id}} = $current;
			$flag = 1;
		}													
		else
		{			
			$current->{sequence} .= uc"$line";
		}	
	}
	close INPUT;
	my $t_finish_get_query_fasta = new Benchmark;

# get and store query sequence lengths  
	foreach my $next (keys %seqs)
	{
		$seqs{$next}->{seq_length} = length ($seqs{$next}->{sequence});
	}

 # open blast file, check each line, keep only those above threshold
 # fix any id_numbers for subject matches that may have been truncated by fastacmd
 	my $pct_query_coverage; 	
 	my $pct_match_coverage;
 	my @suffixlist = ("faa", "fasta", "fsa", "m8", "blast", "filt60");
 	my $basename = fileparse($blast_tab, @suffixlist);
 	my $extension = $pct_cutoff * 100;
 	my $blast_outfile = "$basename"."filt"."$extension";
 	my $query_outfile = "$blast_outfile"."_query_info";
 	my $match_outfile = "$blast_outfile"."_match_info";
 	my %uniq_query_list = ();
	my %uniq_match_list = ();

 	open (BLAST, ">$blast_outfile") or die "can't open filtered blast result $blast_outfile\n $!\n";
 	open (QINFO, ">$query_outfile") or die "can't open query info_file $query_outfile\n $!\n";
 	open (MINFO, ">$match_outfile") or die "can't open match info_file $match_outfile\n $!\n";
 	open (INPUT,"$blast_tab") or die "can't open $blast_tab\n$!\n";	
 	while (<INPUT>) 
	{
	# skip blank lines and header	
		next if ($_ =~ /^\s/);
		next if $_ =~ /query_id/;
		next if $_ =~ /query/i;
		chomp;
	# get rid of extra spaces inserted by ncbi blast
		$_ =~ s/ //g;
		
	# parse entries	
		my @tmp = split "\t", $_;
		my ($query_id, $subject_id, $pct_id, $align_length, $mismatches, $gap_opens, $q_start, $q_end, $s_start, $s_end, $evalue, $bit_score)
			= split "\t", $_;
	
 		my $match_info = $match_lengths{$subject_id};
 		unless (defined $match_lengths{$subject_id})
 		{
 			print STDERR "  WARNING: Blast match subject_id not found in darhorse database\n  $_\n"; 
 		}
 		my ($ncbi_tax_id, $species_name, $match_len, $full_db_id) = split "\t", $match_info;	 		 		
		my $query_len = $seqs{$query_id}->{seq_length};

		unless (defined $seqs{$query_id}->{seq_length} &&  $seqs{$query_id}->{seq_length} >0)
		{
			die "can't find seq_length for query $query_id\n"; 
		}	
		unless ( $match_len >0)
		{
			warn "  WARNING: can't find seq_length for match $subject_id ($full_db_id)\n";
			$match_len = 1;
		}
			
		$pct_match_coverage = $align_length /$match_len;			
		$pct_query_coverage = $align_length / $query_len;
		
	# correct subject_id to full_id corresponding to gi_num so it matches db_seq_id database table	
	
	# remove blast matches that don't meet alignment length criteria	
		if ($pct_query_coverage > $pct_cutoff && $pct_match_coverage > $pct_cutoff)
		{
			my $outstring = join "\t", ($query_id, $full_db_id, $pct_id, $align_length, 
				$mismatches, $gap_opens, $q_start, $q_end, $s_start, $s_end, $evalue, $bit_score);
			print BLAST "$outstring\n";			
			
	# match info file includes info for match sequence only once in, even if hits multiple queries
	# substitute ! for pipe character temporarily to work around perl hash exists bug
			$full_db_id =~ s/\|/!/g;
			unless (exists $uniq_match_list{$full_db_id})
			{
				$uniq_match_list{$full_db_id}++;
				$full_db_id =~ s/!/\|/g;
				my $outstring = join "\t", ($full_db_id, $ncbi_tax_id, $species_name, $match_len);
				print MINFO "$outstring\n";
			}
			
			$query_id =~ s/\|/!/g;
			unless (exists $uniq_query_list{$query_id})
			{
				$uniq_query_list{$query_id}++;
				$query_id =~ s/!/\|/g;			
				print QINFO "$query_id\t$seqs{$query_id}->{seq_length}\t$seqs{$query_id}->{annotation}\n";				
			}			
		}
	}
			
	close INPUT;
	close BLAST;
	close QINFO;
	close MINFO;
	$dbh->disconnect();
	print STDERR "Blast filtering complete\n";
	
	
###################################################
# SUBROUTINES
###################################################
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

sub check_args
{

	if($USAGE || !$query_fasta || !$blast_tab || !$config_filename) 
    {
		print $message;
		exit(0);
	}
	unless (-s $blast_tab)
	{
		print "cannot open blast INPUT $blast_tab\n";		
		exit (0);
	}	
	
	unless (-s $query_fasta)
	{
		print "cannot open query fasta INPUT $query_fasta\n";		
		exit (0);
	}
	
	unless ($pct_cutoff >0.000000001 && $pct_cutoff < 1.001)
	{
		die "Illegal value for filter_threshold $pct_cutoff (must be between 0 and 1)\n";
	}
	unless ($max_lines_per_packet =~ /^\d+$/)
	{
		print STDERR "Illegal value for max_lines_per_packet:$max_lines_per_packet\n";
		exit(0);	
	}	
}

sub new_seq {
  my ($className, $param) = @_;
  my $self = {};
  bless $self, $className;
  my @properties = @$param;
  my $header = join " ", @properties;
  $self->{header} = $header;
  
  $properties[0] =~ s/\>//;
  $self->{id} = $properties[0];
  shift @properties;
  $self->{annotation} = join " ", @properties;
  
  unless (defined $self->{id})
  {
  	warn "no id found for $self->{header}\n";
  }
  
  $self->{seq_length} ="";
  
  return($self)
}

sub sql_select
{		
	my ($list, $table_name, $firstline) = @_;
	my @select_list = @$list;
	my $total_num = scalar @select_list;
	my $sql = $firstline;
	my @row = ();
	my $element = ();
	
	my @results_tbl = ();
		
	unless ($total_num >0)
	{
		print STDERR "couldn't find any selections in list for $table_name\n";
		exit(0);
	}
	my $count = 0;			
	foreach (@select_list)
	{
		$count++;
		chomp;
		$_ =~ s/'//g;		# get rid of internal quotes
		$sql .= "'$_',";
	# send large select statements in pieces to avoid exceeding MySql max_allowed_packet size	
		if ($count % $max_lines_per_packet == 0)
		{
			$sql = substr($sql, 0, -1); # get rid of extra newline and comma, if present
			$sql .= ");\n";
						
			$sth = $dbh->prepare($sql)
    			or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
			$sth->execute()
    			or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");
			@row = ();
			while (@row = $sth->fetchrow_array)
			{				
				my $row_txt;
				foreach my $element (@row)
				{
					$row_txt .= "$element\t";
				}
				push @results_tbl, "$row_txt";
			}
			$sql = $firstline;
		}
		elsif ($count == $total_num)
		{
			$sql = substr($sql, 0, -1); # get rid of extra newline and comma, if present
			$sql .= ");\n";
			$sth = $dbh->prepare($sql)
    			or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
			$sth->execute()
    			or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");
			@row = ();
			while (@row = $sth->fetchrow_array)
			{
		 		my $row_txt;
				foreach my $element (@row)
				{
					$row_txt .= "$element\t";
				}
				push @results_tbl, "$row_txt";
			}
		}
	}	
	return @results_tbl;
}
