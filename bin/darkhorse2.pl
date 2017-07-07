#!/usr/bin/env perl
# darkhorse2.pl 
# Sheila Podell
# May 21, 2017 

##############################################
# PURPOSE

# Script to calculate phylogenetic relatedness of blastp hits to a set of query sequences
# Quantifies blastp match relatedness with a lineage probability score

# General Algorithm:

# Get blastp search against nr (or custom) database in tab-delimited format
# Exclude hits against "self" (user supplied exclusion list of name(s) and/or tax ids) 
# Build a set of guide probabilities based on lineage of tophits
# Use the guide probabilities to identify "corrected" equivalent tophits
# outputs tab-delimited file showing match properties

################################################
# INPUT DATA 

#	tab-delimited blast result file 
	# diamond BLASTP tab | ncbi blast+ -outformat 6 | old ncbi blast blastall m8 format
#	list of one or species names and/or tax ids to exclude
# 	threshold values (bitscore comparison parameter for determining match "equivalency")
#	query (amino acid) fasta file used for blast search (to get annotations and sequence lengths)
#	accessory databases implemented in MySQL, linking reference sequence ids and taxonomy
#   optional: query_info and match_info files
#		shortens time for runs against same blast data set

################################################
# ACCESSORY SCRIPTS REQUIRED

#   get_seq_info.pl
#   getseq_multiple.pl
#   select_tab_first.pl
#   select_tab_equivalent.pl
#   calc_lineage_probability.pl
#   filter_blast.pl

################################################
# OUTPUT FILES

# This program stores initial calculations in files, rather than keeping
# in local variables, to avoid crashing due to memory constraints when scaled
# to full size - routine input data will be hundreds of thousands to 
# millions of lines. These are erased unless debug option is enabled.

# Output includes a logfile (containing warnings and errors), a tab-delimited
# overall summary file, and tally files by species and genus

################################################
use warnings; 
use strict;
use Getopt::Long;
use Cwd;
use File::Basename;
use File::Spec;
use DBI;
use Benchmark;		# for timing code 

# Set default values
	my $t0 = new Benchmark;
	my $tab_input; # blast data in ncbi m8 format";
	my $exclude_listfile; # = "exclude file name;
	my $filter_threshold = 0.1;	
	my $logfilename = "$$"."_logfile";
	my ($program_directory, $script_path, $db_name, $db_program, $db_host, $db_user, 
		$db_user_password, $max_lines_per_packet, $min_lineage_terms, $min_aln_coverage,
		$config_filename, $query_info_file, $match_info_file, $blast_filter);
	my %query_list = (); # key = query_ID, value = object_ref 
	my $debug = 0; # 0 = off 1 = on
	my $genome_fasta = "";
	my ($table_name, $firstline, $getlist);
	my ($sql, $sth, $select_names);	
	my @getlist = ();
	my @row = ();
	my $element = "";	
	
# get command line arguments, check validity
	my $USAGE;
	my $message = qq(
  Usage: $0 -c config file -t blast_tab_infile -e exclude list -g self.fasta 
	
  Required parameters
  	-t  <tab delimited blast input file>
  	-e  <exclude list>
  	-g  <query sequence source file (fasta format)>
  	-c  <config path/file name>
 
  Optional parameters:
	-f  <filter threshold> [default = 0.1] 
	-n  <minimum number lineage terms>
	-d  <debug>
	-b  <blast filter>  [default = 0.7]
	-q  <query_info_file>
	-m  m<atch_info_file>
	
  Description:
	Defines phylogenetic relatedness of blastp hits for a
	set of proteins against NCBI Genbank nr database 
	using a lineage probability score
	)."\n";

	GetOptions( "t=s" => \$tab_input,		# blast result path/filename					
				"e=s" => \$exclude_listfile,	# list of excluded species names and/or tax id's
				"g=s" => \$genome_fasta,	# path/filename for query fasta sequences
				"f=f" => \$filter_threshold,	# filter threshold for equivalence
        		"c=s" => \$config_filename,
        		"n=i" => \$min_lineage_terms,
        		"d=i" => \$debug,
        		"b=f" => \$blast_filter,           #  minimum % coverage for blast alignments
        		"q=s" => \$query_info_file,
        		"m=s" => \$match_info_file,        		
        );
        
 # get global parameters from config file
	unless (defined $config_filename)
		{
			print STDERR "No config filename defined\n";
			print STDERR $message;	
			exit(0);	
		}
		
	my %dh_config = &read_config_file($config_filename);
	$program_directory = $dh_config{program_directory}|| "not found";
	$script_path = "$program_directory/bin/accessory_scripts";
	$db_name = $dh_config{db_name} || "not found";
	$db_program = $dh_config{db_program} || "not found";
	$db_host = $dh_config{db_host} || "not found";
	$db_user = $dh_config{db_user} || "not found";
	$db_user_password = $dh_config{db_user_password} || "not found";
	$max_lines_per_packet = $dh_config{max_lines_per_packet} || "not found";
	$min_aln_coverage = $dh_config{min_align_coverage};	
	unless (defined $min_lineage_terms)
		{
			$min_lineage_terms = $dh_config{min_lineage_terms} || "not found";
		} 
	
	&check_args();
	my $result = &check_db_connection(\%dh_config);
    unless ($result ==1)
    {
        print "\n  ERROR - unable to access database specified in configuration file.\n";
        exit(0);
    }
    print "\nDatabase connection OK\n"; 
    
# blast filter if necessary
   my $t_start_blast_filter = new Benchmark;
	# command line definition of blast filter should over-ride config file, if valid
		unless (defined $blast_filter && $blast_filter >=0 && $blast_filter <= 1)
		{
			$blast_filter = $min_aln_coverage;				
		}
   # if config file had zero or bad value, skip filtering
		unless (defined $blast_filter && $blast_filter >=0 && $blast_filter <= 1)
		{
			$blast_filter = 0;
		}    
    #if (defined $blast_filter && $blast_filter >=0 && $blast_filter <= 1 )
    if ($blast_filter > 0)
    {
        print STDERR "\nFiltering blast input at $blast_filter alignment coverage\n";
		my $pct_cutoff =$blast_filter;
		`$script_path/filter_blast.pl -q $genome_fasta -t $tab_input -c $config_filename -f $pct_cutoff`;
  
		my @suffixlist = ("faa", "fasta", "fsa", "m8", "blast", "blastp");
		my $basename = fileparse($tab_input, @suffixlist);
		my $extension = $pct_cutoff * 100;
		my $blast_outfile = "$basename"."filt"."$extension";
		my $query_outfile = "$blast_outfile"."_query_info";
		my $match_outfile = "$blast_outfile"."_match_info";
		$tab_input = $blast_outfile;
	
		unless (defined $query_info_file && -s $query_info_file)
		{
			$query_info_file = $query_outfile;
		}
	
		unless (defined $match_info_file && -s $match_info_file)
		{
			$match_info_file = $match_outfile;
		}   
    }
    else
    {
    	print STDERR "\n  WARNING: Blast filter setting = $blast_filter. No filtering performed for minimum alignment coverage\n";
    }
    my $t_finish_blast_filter = new Benchmark;
    
# convert relative paths to absolute 
	my ($tab_file, $tab_dir) = fileparse($tab_input);
	my ($fasta_file, $fasta_dir) = fileparse($genome_fasta);
	my ($exclude_file, $exclude_dir) = fileparse($exclude_listfile);
	my ($program_file, $program_dir) = fileparse($0);	
	$tab_dir = File::Spec->rel2abs($tab_dir);
	$fasta_dir = File::Spec->rel2abs($fasta_dir);
	$exclude_dir = File::Spec->rel2abs($exclude_dir);	
	$tab_input = "$tab_dir/$tab_file";
	$genome_fasta = "$fasta_dir/$fasta_file";
	$exclude_listfile = "$exclude_dir/$exclude_file";
	
	if (defined $match_info_file && -s $match_info_file)
	{
		my ($minfo_file, $minfo_dir) = fileparse($match_info_file);
		$minfo_dir = File::Spec->rel2abs($minfo_dir);
		$match_info_file = "$minfo_dir/$minfo_file";
	}
	
	if (defined $query_info_file && -s $query_info_file)
	{
		my ($qinfo_file, $qinfo_dir) = fileparse($query_info_file);	
		$qinfo_dir = File::Spec->rel2abs($qinfo_dir);
		$query_info_file = "$qinfo_dir/$qinfo_file";
	}

# set up a working directory, move there
	my $result_dir = "calcs"."_"."$$"."_filt_$filter_threshold";
	mkdir $result_dir, 0755 or warn "Cannot make results directory $result_dir\n: $!";
	chdir $result_dir;
	my $cwd = getcwd();

# open a log file
	my $date = `date`;
	chomp $date;
	open (LOGFILE, ">$logfilename") or die "can't open log file $logfilename\n $!\n";
	print LOGFILE qq($date
$program_directory/$program_file
	
Input parameters:
	tab_input_file		$tab_file				
	exclude_file_name 	$exclude_file
	fasta_file_name		$fasta_file
	filter_threshold	$filter_threshold
	config_file			$config_filename\n);	

	if ($debug ==1)
	{
		foreach my $key (sort keys %dh_config)
		{
			print LOGFILE "\t\t$key=$dh_config{$key}\n";
		}
	}
	
# connect to database
	my $db_path = "DBI:$db_program:$db_name:$db_host:";
	my $dbh = DBI->connect($db_path, $db_user, $db_user_password) or die "Cannot connect: $DBI::errstr";

##############################
# Even if sequences with fewer than minimum number of terms have been pre-excluded from
# reference fasta file used for initial blastp search, need to deal with case where
# database creation parameters are less stringent than minumum number of terms in
# current analysis

# Get list of database id numbers with missing lineage information
	my $t_start_missing_lineage_search = new Benchmark;
	print STDERR "\nFinding database entries with insufficient lineage information (<$min_lineage_terms terms)\n";
	$table_name = "lineage_index";
	
	my @bad_ids = ();
	$sql = "SELECT full_id FROM lineage_index where num_terms < $min_lineage_terms;";
	$sth = $dbh->prepare($sql)
    	or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
	$sth->execute()
		or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");
	while (@row = $sth->fetchrow_array)
	{				
		my $row_txt;
		foreach my $element (@row)
		{
			$row_txt .= "$element\t";
		}
		push @bad_ids, "$row_txt";
	}
	
# put sql results into a hash
	my %missing_lineage_ids = ();
	my $missing_lineage_count = 0; 

	foreach (@bad_ids)
	{
		chomp;
		my $full_id = $_;
	# get rid of extra spaces and tabs inserted by blast/diamond software
		$full_id =~ s/ //sg;
		$full_id =~ s/\t//sg;
		$missing_lineage_ids{"$full_id"} = 1;
		$missing_lineage_count++;
	}	

	my $t_finish_missing_lineage_search = new Benchmark;
	
	if ($missing_lineage_count < 1)
	{
		#print STDERR "Warning: didn't find any id numbers for missing lineages\n";
	}

	if ($debug ==1)
	{
		print LOGFILE "\n\nFound $missing_lineage_count database entries without lineage information\n";
	}
	
# Get species name, tax_id and lineage for all matches from SQL database 
	print STDERR "Finding species names for blast matches\n";
	my $t_start_find_species = new Benchmark;
	$getlist = `cut -f2 $tab_input`; 
	@getlist = split "\n", $getlist;
	
# remove any match id numbers that have no lineage data from list
# also remove any matches that have less than the minumum number of terms for species names
 	my @filtered_getlist = ();
	foreach (@getlist)
 	{		
		chomp;		
		if (exists $missing_lineage_ids{$_})
 		{
 			next;
 		}
		push @filtered_getlist, $_;
 	}

	my $num_removed = (scalar @getlist) - (scalar @filtered_getlist);
	if ($debug ==1)
	{
		print LOGFILE "Found $num_removed matches with insufficient lineage information\n\n";
	}
	
	$table_name = "db_seq_ids";
	$firstline = "SELECT full_id, ncbi_tax_id, species_name FROM $table_name where full_id in (";	
	
# check to see if match_info species table already exists. If not, create it
	my @species_table = ();
	 if (defined $match_info_file && -s $match_info_file)
 	{
 		print LOGFILE "\tmatch_info source \t$match_info_file\n";
 		open (INPUT, "<$match_info_file") or die "can't open $match_info_file\n $!\n";
 		while (<INPUT>)
 		{
 			push @species_table, $_;
		}
 		close INPUT;		
 	}
 	else
	{
		print STDERR "Getting match info from database\n";
		@species_table=&sql_select(\@filtered_getlist,$table_name, $firstline);
	}
	
# get lineage lookup table for all tax_ids found among matches
	my %tax_ids = ();
	foreach my $entry (@species_table)
	{
		my @tmp = split '\t', $entry;
		my $lookup_id = $tmp[1];
		next unless $lookup_id =~ /^-*\d+$/; # accomodate negative taxid numbers as well as positive
		$tax_ids{$lookup_id}++;
	}
	my @tax_id_list1 = keys %tax_ids;
	$table_name = "precalc_lineages";
	$firstline = "SELECT ncbi_tax_id, lineage from $table_name where ncbi_tax_id in (";
	
	my @lineage_lookup=&sql_select(\@tax_id_list1,$table_name, $firstline);
	foreach (@lineage_lookup)
	{
		my @tmp = split '\t', $_;
		$tax_ids{$tmp[0]}=$tmp[1];
	}	
	
	my $t_finish_find_species = new Benchmark;

# Find match id numbers corresponding to self, to be excluded
	my $t_start_find_self = new Benchmark;	
	print STDERR "Removing hits to self-id numbers.\n";
	print LOGFILE "\nExcluding the following terms:\n";
	
	my %exclude_terms = &get_exclude_terms($exclude_listfile);
	foreach (sort keys %exclude_terms)
	{
		print LOGFILE "\t$_\t($exclude_terms{$_})\n";
	}
	my %self_ids = &self_filter(\%exclude_terms, \@species_table);	

	my $self_count = scalar (keys %self_ids);
 	
	my $t_finish_find_self = new Benchmark;

# Process master file: remove matches to self or missing lineages 		
	my $t_start_self_filter = new Benchmark;	

	my $self_exclusion_count = 0;
	$missing_lineage_count = 0;
	my %query_hit_tally = ();
	
	open (INPUT, "<$tab_input") or die "can't open $tab_file\n $!\n";
	open (OUTPUT, ">master.tab.1") or die "can't open outfile master.tab.1\n $!\n";
	
	while (<INPUT>)
	{
		chomp;
		next if ($_ =~ /^\s*$/);
	# get rid of extra spaces inserted by ncbi blast
		$_ =~ s/ //g;
		
		my @tmp = split "\t", $_;
		my $query_id = $tmp[0];
		my $match_id = $tmp[1];
		next unless (defined $match_id && length $match_id >1);						
		unless (exists $self_ids{$match_id}
			|| exists $missing_lineage_ids{$match_id})	
		{
			print OUTPUT "$_\n";
			$query_hit_tally{$query_id}++;
		}		
	}		
	close INPUT;
	close OUTPUT;
	
	#print STDERR "Removed $self_exclusion_count self matches and $missing_lineage_count missing lineage matches.\n";
 	print LOGFILE "\nRemoved $self_exclusion_count self matches and $missing_lineage_count missing lineage matches.\n";

 	my $t_finish_self_filter = new Benchmark;

# Create object to hold data for each query that has non-self hits
	my $t_start_query_info = new Benchmark;
	my @nonself_tally = `cut -f1 master.tab.1  | sort | uniq -c`;
	foreach my $query (@nonself_tally)
	{
		chomp $query;
		$query =~ /\s*(\d+)\s+(.+)/;
		my $query_hits = $1;
		my $query_id = $2;
		
		unless (defined $query_hits && defined $query_id)
		{
			print STDERR "problem with query $query in nonself tally\n";
			print LOGFILE "problem with query $query in nonself tally\n";
			next;
		}
		my @tmp = ($query_id, "not_analyzed", $query_hits);
		my $current = &new_query("query", \@tmp);
		$query_list{$current->{query_id}}= $current;
	}
	
# get sequence length and annotation information for each query object	
	my  @query_hit_info = ();
	if (defined $query_info_file && -s $query_info_file)
	{
		print LOGFILE "query_info read from $query_info_file\n";
		open (INPUT, "<$query_info_file") or die "can't open $query_info_file\n $!\n";
		while (<INPUT>)
		{
			chomp;
			next if ($_ =~ /^\s*$/);
			my @fields = split /\t/, $_;
			next unless (scalar @fields ==3 || scalar @fields ==2);	#allow for queries with no annotation	
			push @query_hit_info, $_;
		}
		close INPUT;
	}	
	else
	{
		print STDERR "Getting query info from fasta file\n";
		`cut -f1 master.tab.1 |sort| uniq > query_ids`;
		`$script_path/getseq_multiple.pl query_ids $genome_fasta > query_hits.fasta`;
		@query_hit_info = `$script_path/get_seq_info.pl query_hits.fasta`;	
	}
	my $difference = abs(scalar @query_hit_info - scalar keys %query_list);
	unless ($difference == 0)
	{
		if ($debug ==1)
		{
			print STDERR "Warning: couldn't find detailed query information for $difference queries\n";
			print LOGFILE "Warning: couldn't find detailed query information for $difference queries\n";
		}
	}
	
	foreach my $line (@query_hit_info)
	{
		chomp $line;
		next if ($line =~ /^\s*$/);
		my @tmp = split "\t", $line;
		my $query_id = $tmp[0];
		
	# make sure this query had valid matches before adding object info to query list 	
		if (exists $query_list{$query_id})
		{
			$query_list{$query_id}->{seq_length} = $tmp[1];
			$query_list{$query_id}->{annotation} = $tmp[2] || "no data";
			$query_list{$query_id}->{nonself_tally} = $query_hit_tally{$query_id} || "no data";
		}
	}
	
	# How many queries had non-self mataches?
		my $num_genome_proteins = `grep -c '\>' $genome_fasta`;
		chomp $num_genome_proteins;
		
		my $num_matched_proteins = scalar (keys %query_list);#@query_hit_info;
		
		print STDERR  "$num_matched_proteins/$num_genome_proteins query proteins had non-self matches.\n";
		print LOGFILE "$num_matched_proteins/$num_genome_proteins query proteins had non-self matches.\n";
	
	unless ($debug > 0)
	{
		unlink "query_ids";
		unlink "query_hits.fasta";
	}
	my $t_finish_query_info = new Benchmark;

# Filter master file tab output by bitscores (select equivalent)
	my $t_start_threshold_filter = new Benchmark;
	my $percent = $filter_threshold * 100;
	print STDERR "Selecting \"equivalent\" hits (relative bitscore within $percent%).\n";
		`$script_path/select_tab_equivalent.pl -i master.tab.1 -k 0 -s 11 -p pct -m max -t $filter_threshold -l 500 > master.tab.2`;
		
	my $filt1 = &count_lines("master.tab.1");
	my $filt2 = &count_lines("master.tab.2");
	print LOGFILE "Equivalence filter reduced blast hits from $filt1 to $filt2 \"besthits\".\n";
	
	unless ($debug > 0)
	{
		unlink "master.tab.1";
	}

# count number of candidates for each query
	my @candidates_tally = `cut -f1 master.tab.2  | sort | uniq -c`;
	my $max_set_size = 0;
	my $conserved_query = "";
	foreach my $query (@candidates_tally)
	{
		chomp $query;
		$query =~ /\s+(\d+)\s+(.+)/;
		my $query_hits = $1;
		my $query_id = $2;
		if ($query_hits >$max_set_size)
		{
			$max_set_size = $query_hits;
			$conserved_query = $query_id;
		}
		
		unless (defined $query_hits && defined $query_id)
		{
			print STDERR "problem with query $query in nonself tally\n";
			print LOGFILE "problem with query $query in nonself tally\n";
			next;
		}
		$query_list{$query_id}->{candidates_tally} = $query_hits;
	}
	my $short_annot = substr ($query_list{$conserved_query}->{annotation}, 0, 60);
	print LOGFILE "Maximum candidate set size was $max_set_size\n";
	print LOGFILE "Most conserved query: $conserved_query $short_annot\n";
	
	my $t_finish_threshold_filter = new Benchmark;
	
# Select (first round) tophits for each query to get guide probabilities
	print STDERR "Getting first round LPI guide scores\n";
	my $t_start_get_lineages = new Benchmark;
	# make sure the tab file is sorted by query name, then descending bitscore 
		` sort -k 1,1 -k 12,12nr master.tab.2 > master.tab.3`;
		my @first_round_tophit_lines = `$script_path/select_tab_first.pl master.tab.3`;
		unless ($debug > 0)
		{
			unlink "master.tab.2";
		}
	 		
	# get match id numbers for tophits, add info to query object
		my %first_round_lineages = ();
		foreach (@first_round_tophit_lines)
		{
			next if ($_ =~ /^\s+$/);
			my @tmp = split "\t", $_;
			my $query_id = $tmp[0];
			my $match_id = $tmp[1];
			$query_list{$query_id}->{tophit_line} = $_;
			$query_list{$query_id}->{tophit_id} =$match_id;
			$first_round_lineages{$match_id} = "";					
		}
		
	# get species names and lineages for each tophit
		@getlist = (keys %first_round_lineages);
		my %guide_probabilities = &calc_LPI(\@getlist); #key = lineage #value = LPI
		
if ($debug > 0)
{
	open DEBUG, ">guide_probabilities";
	foreach my $key (keys %guide_probabilities)
	{
		print DEBUG "$key\t$guide_probabilities{$key}\n";
	}
	close DEBUG;			
}
# Append species names, lineages, guide probabilities to master table
# For species that didn't make it into first round tophits, guide probability = 0
	my $matchnames = `cut -f2 master.tab.3`;
	@getlist = split "\n", $matchnames;
	$table_name = "lineage_index";
	$firstline = "SELECT full_id, species_name, ncbi_tax_id, lineage FROM $table_name where full_id in (";
	my @filtered_match_info = &sql_select(\@getlist, $table_name, $firstline);
	my %filtered_match_addlines =();
	foreach (@filtered_match_info)
	{
		my @tmp = split "\t", $_;
		my $full_id = shift @tmp;
		$filtered_match_addlines{$full_id} = join "\t", @tmp;
		#my $length  = length $filtered_match_addlines{$full_id};
	}	
	
	open (INPUT, "<master.tab.3") or die "can't open master.tab.3\n $!\n";
	open (OUTPUT, ">master.tab.4") or die "can't open outfile master.tab.4\n $!\n";

	while (<INPUT>)
	{
		chomp;
		my @tmp = split "\t", $_;
		my $match_id = $tmp[1];
		
		unless (exists $filtered_match_addlines{$match_id})
		{
			my $line_num = $.;
			my $msg = "\n  BLAST INPUT ERROR.\n";
			$msg .= "  Couldn't find match id \"$match_id\" in database $db_name\n";
			$msg .=  "$_\n"; 
			print STDERR "$msg";
			print LOGFILE "$msg";
			next;	
		}
				
		unless (length $filtered_match_addlines{$match_id} > 0)
		{
			print STDERR "no filtered_match_addlines for \n$_\n";
		}
		
		my $modified_line = "$_"."\t"."$filtered_match_addlines{$match_id}";
		@tmp = split "\t", $modified_line;
		my $lineage = $tmp[-1];
		
		# debug
		if ($debug >0)
		{
			unless (defined $lineage && length $lineage >0)
			{
				print "undefined lineage line 665\t$modified_line\n";
			}
		}
			
		unless (exists $guide_probabilities{$lineage})
		{
			$guide_probabilities{$lineage} = 0;	
		}
		$modified_line .= "\t$guide_probabilities{$lineage}";
		print OUTPUT "$modified_line\n";
	}
	close INPUT;
	close OUTPUT;
	unless ($debug > 0)
	{
		unlink "master.tab.3";
	}

# Select (second round) tophits for each query 
	print STDERR "Getting second round (final) LPI scores\n";
	# Sort master file by query ID, then decreasing lineage probability (field 16),
	# then by decreasing bitscore (field 12)
		`sort -k 1,1 -k 16,16nr -k 12,12nr master.tab.4 >  master.tab.5`;
		my @second_round_tophit_lines = `$script_path/select_tab_first.pl master.tab.5`;
		unless ($debug > 0)
		{
			unlink "master.tab.4";
		}
		
	# get match id numbers for 2nd round tophits, add info to query object
		my %second_round_lineages = ();
		foreach (@second_round_tophit_lines)
		{
			next if ($_ =~ /^\s+$/);
			my @tmp = split "\t", $_;
			my $query_id = $tmp[0];
			my $match_id = $tmp[1];
			$second_round_lineages{$match_id} = "";
			$query_list{$query_id}->{besthit_id} =$match_id;
			
		# blast data should still be in original m8 order
			$query_list{$query_id}->{besthit_pct_id} = $tmp[2];
			$query_list{$query_id}->{besthit_align_len} = $tmp[3];
			$query_list{$query_id}->{besthit_evalue} = $tmp[10];
			$query_list{$query_id}->{besthit_bitscore} = $tmp[11];
			$query_list{$query_id}->{besthit_species} = $tmp[12];
			$query_list{$query_id}->{besthit_tax_id} = $tmp[13];
			$query_list{$query_id}->{besthit_lineage} = $tmp[14];	
		}
		
		@getlist = (keys %second_round_lineages);
		my %LPI_results = &calc_LPI(\@getlist); #key = lineage #value = LPI	
		my @sorted_lpi_values = sort {$b <=> $a} values %LPI_results;
		my $lpi_max = $sorted_lpi_values[0];		
		my $t_finish_get_lineages = new Benchmark;
		
if ($debug > 0)
{
	open DEBUG, ">second_round_lineages";
	foreach my $key (keys %LPI_results)
	{
		print DEBUG "$key\t$LPI_results{$key}\n";
	}
	close DEBUG;			
}

# get annotations for best hits
	print STDERR "Getting annotations for besthits\n";
	my $t_start_get_besthit_annotations = new Benchmark;
	
	@getlist = (keys %second_round_lineages);
	$table_name = "db_seq_ids";
	$firstline = "SELECT full_id, annotation FROM db_seq_ids WHERE full_id in (";
	my @annotation_output = &sql_select(\@getlist, $table_name, $firstline);	
	
	my %besthit_annotations =();
	foreach (@annotation_output)
	{
		my @tmp = split "\t", $_;
		my $match_id = $tmp[0];
		my $match_annotation = $tmp[1];
		$besthit_annotations{$match_id} = $match_annotation;
	}
	
	my $t_finish_get_besthit_annotations = new Benchmark;
	
# add final best hit information to query objects
	foreach (keys %query_list)
	{
		my $next = $query_list{$_};
		my $besthit_id = $next->{besthit_id};
		my $besthit_lineage = $next->{besthit_lineage};
		
		$next->{besthit_annotation} = $besthit_annotations{$besthit_id};
		$next->{besthit_LPI} = $LPI_results{$besthit_lineage} || "not found";
	}

# print out summary, sorted in order of decreasing LPI, lineage, species
	my $smry_name = "$$"."_smry";
	open (SMRY, ">$smry_name") or die "can't open summary file $smry_name\n $!\n";
	my $headers = "query_id\tnonself_matches\tcandidate_matches\tbesthit_id".
				"\traw_LPI\tnorm_LPI\tpct_id\tquery_len\talign_len".
				"\tpct_coverage\tevalue\tbitscore".
				"\ttax_id\tspecies\tlineage\tquery_annotation".
				"\tbesthit_annotation\n";
	print SMRY "$headers";
		
	foreach (sort {
			$query_list{$a}->{besthit_LPI} 
			<=> $query_list{$b}->{besthit_LPI} 
			
			|| $query_list{$a}->{besthit_lineage} 
			cmp $query_list{$b}->{besthit_lineage}
			
			|| $query_list{$a}->{besthit_species} 
			cmp $query_list{$b}->{besthit_species}
			
			} 
			keys %query_list)
	{	
		my $outstring = &query_toString($query_list{$_}, $lpi_max);
		print SMRY "$outstring\n";	
	}	
	close SMRY;
	
# check for errors
	&check_errs($smry_name);
	
# get summary statistics
	my $t_start_get_stats = new Benchmark;
	
	# species tally
	`cut -f 14 $smry_name | sort | uniq -c | sort -nr > species_tally`;
	
		my %genus_tally = (); # key = genus name # value = number of hits
		open (SMRY, "<$smry_name") or die "can't open summary file $smry_name\n $!\n";
		while (<SMRY>)
		{
			my @fields = split "\t",$_;
			my $species_name = $fields[13];
			my @terms = split " ", $species_name;
			my $genus = $terms[0];			
			if ($species_name =~/Candidatus/ )
			{				
				$genus = "$terms[0]". " "."$terms[1]"; 
			}
			if ($species_name =~/^bacterium/ 
			|| $species_name =~/^uncultured/
			|| $species_name =~/^marine/
			|| scalar @terms < 2)
			{				
				$genus = "$species_name"; 
			}			
			$genus_tally{$genus}++;				
		}	
		close SMRY;
		open (GENUS_TALLY, ">genus_tally") or die "can't open summary file genus_tally\n $!\n";
		foreach my $key (sort {$genus_tally{$b} <=> $genus_tally{$a}} keys %genus_tally)
		{
			print GENUS_TALLY "$genus_tally{$key}"."\t"."$key\n";
		}
		close GENUS_TALLY;
	
sub new_genus {
  my ($className, $param) = @_;
  my $self = {};
  bless $self, $className;
  my @properties = @$param;
   
  $self->{full_name} = $properties[0];
  $self->{seq_length} = 0;
  $self->{num_blasthits} = $properties[1] || "no data";
  $self->{nonself_tally} = $properties[2] || "no data";;
  $self->{candidates_tally} = 0;
		   
  return($self)
}
# LPI histogram
	my $histogram_min = 0;
	my $histogram_max = 1.0;
	my $num_bins = 50;	
	
# put raw LPI histogram into file
	my @lpi_lines = `cut -f5 $smry_name`;
	my @lpi_scores = ();
	foreach my $line (@lpi_lines)
	{			
		chomp $line;
		next unless ($line =~ /^[\d\.]+$/);
		push @lpi_scores, $line;
	}
	my %lpi_histogram = &histogram_calcs($histogram_min, $histogram_max, $num_bins, \@lpi_scores);	
	my $lpi_histogram_file = "raw_lpi_histogram";
	open (OUTPUT, ">$lpi_histogram_file") or die "can't open output file $lpi_histogram_file\n $!\n";

	my @sorted_keys = sort keys %lpi_histogram;
	foreach my $key (@sorted_keys)
	{
		print OUTPUT "$key\t$lpi_histogram{$key}\n";
	
	}
	close OUTPUT;
	
# put normalized LPI histogram into file
	my @norm_lpi_lines = `cut -f6 $smry_name`;
	my @norm_lpi_scores = ();
	foreach my $line (@norm_lpi_lines)
	{			
		chomp $line;
		next unless ($line =~ /^[\d\.]+$/);
		push @norm_lpi_scores, $line;
	}
	my %norm_lpi_histogram = &histogram_calcs($histogram_min, $histogram_max, $num_bins, \@norm_lpi_scores);	
	my $norm_lpi_histogram_file = "norm_lpi_histogram";
	open (OUTPUT, ">$norm_lpi_histogram_file") or die "can't open output file $norm_lpi_histogram_file\n $!\n";

	@sorted_keys = sort keys %norm_lpi_histogram;
	foreach my $key (@sorted_keys)
	{
		print OUTPUT "$key\t$norm_lpi_histogram{$key}\n";
	
	}
	close OUTPUT;
	
	my $t_finish_get_stats = new Benchmark;
	
# clean up, Check time benchmarks
	my $t_final = new Benchmark;
	$dbh->disconnect();
	unless ($debug >0)
	{
		unlink "master.tab.5";
	}

	my $t_total = timediff ($t_final, $t0);
	my $t_blast_filter = timediff($t_finish_blast_filter, $t_start_blast_filter);
	my $t_get_species = timediff ($t_finish_find_species, $t_start_find_species);
	my $t_find_self = timediff ($t_finish_find_self, $t_start_find_self);
	my $t_self_filter = timediff ($t_finish_self_filter, $t_start_self_filter);
	my $t_query_info = timediff ($t_finish_query_info, $t_start_query_info);
	my $t_threshold_filter = timediff ($t_finish_threshold_filter, $t_start_threshold_filter);
	my $t_get_lineages = timediff ($t_finish_get_lineages, $t_start_get_lineages);
	my $t_get_besthit_annotations = timediff ($t_finish_get_besthit_annotations, $t_start_get_besthit_annotations);
	my $t_get_stats = timediff ($t_start_get_stats, $t_finish_get_stats);
	
	
	print LOGFILE "\n\nTotal time elapsed: ", timestr($t_total), "\n";
	
	if (defined $blast_filter)
	{
	    print LOGFILE "\tBlast filter: ", timestr($t_blast_filter), "\n";
	}
	print LOGFILE "\tGet match species: ", timestr($t_get_species), "\n";
	print LOGFILE "\tFind self id numbers: ", timestr($t_find_self), "\n";
	print LOGFILE "\tSelf-filter: ", timestr($t_self_filter), "\n";
	print LOGFILE "\tGet query info: ", timestr($t_query_info), "\n";
	print LOGFILE "\tThreshold filter: ", timestr($t_threshold_filter), "\n";
	print LOGFILE "\tCalc LPIs (2 rounds): ", timestr($t_get_lineages), "\n";
	print LOGFILE "\tGet besthit annotations: ", timestr($t_get_besthit_annotations), "\n";
	print LOGFILE "\tGet statistics: ", timestr($t_get_stats), "\n\n";
	
	close LOGFILE;
	
#########################################################
# SUBROUTINES
#########################################################
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
	if($USAGE || !$tab_input || !$exclude_listfile || !$genome_fasta) 
    {
		print STDERR $message;
		exit(0);
	}	
	unless (-s $tab_input)
	{
		print STDERR "cannot open infile $tab_input\n";		
		exit (0);
	}	
		unless (-s $exclude_listfile)
	{
		print STDERR "cannot open infile $exclude_listfile\n";		
		exit (0);
	}	
	unless (-s $genome_fasta)
	{
		print STDERR "cannot open query source fasta file $genome_fasta\n";		
		exit (0);
	}		
	unless ($filter_threshold >=0 && $filter_threshold <= 1)
	{
		print STDERR "invalid filter threshold value $filter_threshold <range 0-1>\n";
		exit (0);
	}
	unless ($debug == 0 || $debug == 1)
	{
		print STDERR "Debug parameter (-d) must be set at zero (off) or one (on)\n";
		exit(0);
	}
	unless ($max_lines_per_packet =~ /^\d+$/)
	{
		print STDERR "Illegal value for max_lines_per_packet:$max_lines_per_packet\n";
		exit(0);	
	}
	unless ($min_lineage_terms =~ /^\d+$/)
	{
		print STDERR "Illegal value for min_lineage_terms:$min_lineage_terms\n";
		exit(0);	
	}
	if (defined $blast_filter)
	{
	    unless ($blast_filter >=0 && $blast_filter <= 1)
	    {
	        print STDERR "Illegal value for blast filter $blast_filter <range 0-1>\n";
	    
	    }
	
	}
	
# check for required accessory scripts
        my @program_list = ("select_tab_first.pl",
                            "select_tab_equivalent.pl",
                            "getseq_multiple.pl",
                            "get_seq_info.pl",
                            "calc_lineage_probability.pl",
                            "filter_blast.pl",);         
         foreach (@program_list)
         {
           my $program =  "$script_path/$_";
            unless (-s $program && -x $program)
            {
                print STDERR "Couldn't find program $_\n";
                print STDERR "script path = $script_path\n";
                print STDERR "program directory = $program_directory\n";
                exit (0);
            }         
         }	
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
sub sql_select
{		
	my ($list, $table_name, $firstline) = @_;
	my @select_list = @$list;
	my $total_num = scalar @select_list;
	my $sql = $firstline;
	
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

sub get_exclude_terms
{
	my ($exclude_filename) = @_;
	my %exclude_terms = ();
	my $exclude_type = "name";
	
	open (INFILE, $exclude_filename)  or die "can't open input file $exclude_filename\n$!\n";
	while (my $line = <INFILE>) 
	{
		next if $line =~ /^\s+$/;
		chomp $line;
				
	# decide if it is a number or a name
		if ($line =~ /^[-\d]+$/)
		{
			$exclude_type = "number";
		}
		else
		{
			$exclude_type = "name";
		}
		$exclude_terms{$line}= $exclude_type;		
	}	
	close INFILE;
	
	return %exclude_terms;
}

sub self_filter
{
	my ($exclude_hash, $table_ref) = @_;
	my %terms = %$exclude_hash;
	my @species_id_tbl = @$table_ref;
	
	#store self-ids in a hash to ensure each unique id is only included once
	my %exclude_ids_uniq = ();				
	my $flag = 0;
	foreach my $tbl_entry (@species_id_tbl)	
	{
		 my @tmp = split "\t", $tbl_entry;
		 my $match_id = $tmp[0];
		 my $tax_id = $tmp[1];
		 my $species_id =$tmp[2];
		 my $lineage = $tax_ids{$tax_id};
		 $tbl_entry .= "\t$lineage";

	# exclude any match sequences with no species information in db 		
		unless (defined $species_id && !($species_id =~ /null/i) )
		{
			$exclude_ids_uniq{$match_id} = "no species excl\n";		
		}		

	# check tax_ids - make sure numbers match EXACTLY. words can match anywhere.
		foreach my $key (keys %terms)
		{						
			if ($tbl_entry =~ /$key/i)
			{
				unless ($terms{$key} eq "number")
				{
					$exclude_ids_uniq{$match_id} = "self name excl";
				}
				if ($terms{$key} eq "number" && $tax_id =~ /^$key$/) # && $tbl_entry =~ /^$key$/)
				{					
					$exclude_ids_uniq{$match_id} = "self number excl";
				}
			}
		}
	}	
	return (%exclude_ids_uniq);
}

sub new_query {
  my ($className, $param) = @_;
  my $self = {};
  bless $self, $className;
  my @properties = @$param;
  
  unless (defined $properties[0] && length $properties[0] > 2)
  {
  	print STDERR "new_query failed - no ID number\n";
  	print STDERR "input params: @properties\n";
  	exit(0);
  }
  
# quantitative params should be initialized with a number
# string params should be initialized with a string

  $self->{query_id} = $properties[0];
  $self->{seq_length} = 0;
  $self->{num_blasthits} = $properties[1] || "no data";
  $self->{nonself_tally} = $properties[2] || "no data";;
  $self->{candidates_tally} = 0;
  	
  $self->{besthit_line} = "";
  
  $self->{annotation} = "no data";
  $self->{besthit_annotation} = "no data";
  	
  	$self->{tophit_id} =""  || "no data";
  	$self->{tophit_line} = "";

	$self->{besthit_id} =""  || "no data";		
	$self->{besthit_pct_id} = 0;
	$self->{besthit_align_len} = 0;
	$self->{besthit_align_coverage} = 0;
	$self->{besthit_evalue} = 10;
	$self->{besthit_bitscore} = 0 ;
	$self->{besthit_tax_id} =""  || "no data";
	$self->{besthit_species}  =""  || "no data";
	$self->{besthit_lineage} =""  || "no data";
	$self->{besthit_LPI} = 0;
		   
  return($self)
}

sub count_lines
{
	my ($filename) = @_;
	my $num_lines = `wc -l $filename`;
	$num_lines =~ /(\d+).+/;
	return $1;
}

sub calc_LPI
{
	my ($listref) = @_;
	my @id_list = @$listref;
	my %LPI_hash = ();
	
	#use sql query to get lineages from lineage_index table
		$table_name = "lineage_index";
		$firstline = "SELECT lineage FROM $table_name where full_id in (";
		my @first_round_tophit_lineages = &sql_select(\@id_list, $table_name, $firstline);	
		
	# calculate probability for each lineage, put into a hash
		open (OUTPUT, ">lpi_lineage_list") or die "can't open outfile lpi_lineage_list\n $!\n";
		foreach (@first_round_tophit_lineages)
		{
			print OUTPUT "$_\n";
		}
 	my $output = `$script_path/calc_lineage_probability.pl lpi_lineage_list`; 	
 	my @output_lines = split "\n", $output;
 	
 if ($debug > 0)
 {
 	`$script_path/calc_lineage_probability.pl lpi_lineage_list > lpi_lineage_list_output`;
 }	
 	
	my %probabilities = (); # key = lineage, value = preliminary LPI
	
	foreach (@output_lines)
	{
		my @tmp = split "\t", $_;
		my $lineage = $tmp[0];
		my $LPI = $tmp[1];
		$probabilities{$lineage} = $LPI;
	}
	unless ($debug >0)
	{
		unlink "lpi_lineage_list";
	}
	return %probabilities;
}

sub query_toString
{
	my ($query, $lpi_max) = @_;
	my @properties = ();
	push @properties, "$query->{query_id}" || "no data";
	push @properties,"$query->{nonself_tally}" || "no data";
  	push @properties,"$query->{candidates_tally}" || "no data"; 
	push @properties,"$query->{besthit_id}" || "no data";

	
	my $formatted_LPI = sprintf ("%.3f", $query->{besthit_LPI});	
	push @properties, $formatted_LPI;
	
	my $norm_LPI = $query->{besthit_LPI}/$lpi_max;
	my $formatted_norm_LPI = sprintf ("%.3f", $norm_LPI);
	push @properties, $formatted_norm_LPI;	
	
	push @properties,"$query->{besthit_pct_id}" || "no data";
	push @properties,"$query->{seq_length}" || "no data";
	push @properties,"$query->{besthit_align_len}" || "no data";
 
	if (defined $query->{besthit_align_len} && defined $query->{seq_length}
		&& $query->{besthit_align_len} =~ /\d+/ && $query->{seq_length}  =~ /\d+/
		&& $query->{seq_length} >0)
	{
		$query->{besthit_align_coverage} = $query->{besthit_align_len}/$query->{seq_length};
	}
	else 
	{
		$query->{besthit_align_coverage} = 0;
	}
	my $pct_coverage = $query->{besthit_align_coverage} * 100;
	my $formatted_coverage = sprintf ("%d", $pct_coverage);
	push @properties, $formatted_coverage; 
	
# workaround for perl misinterpretation of quotes around zero e-value as an empty string
	if (defined $query->{besthit_evalue})
	{
		push @properties, $query->{besthit_evalue};		
	}
	else
	{
		push @properties, "no data";
	}
	my $formatted_bitscore = sprintf ("%.0f", $query->{besthit_bitscore});
	
	push @properties,"$formatted_bitscore" || "no data";
	push @properties,"$query->{besthit_tax_id}" || "no data";
	push @properties,"$query->{besthit_species}" || "no data";
	push @properties,"$query->{besthit_lineage}" || "no data";
	push @properties,"$query->{annotation}" || "no data";
	if (defined $query->{besthit_annotation})
	{
		push @properties, $query->{besthit_annotation};		
	}
	else
	{
		push @properties, "no data";
	}	
		
	my $string = join "\t", @properties;
	return $string;	
}

sub check_errs
{
	my ($smry_file) = @_;
	
	open (INFILE, "$smry_file")  or die "can't open input file $$.smry\n$!\n";
	while (my $line = <INFILE>) 
	{
	# skip header line
		next if ($line =~ /query/);
	# no self-hits found
		my @fields = split "\t", $line;
		if ($fields[2] < 1)
		{
			print STDERR "Warning: no self hits found for query $fields[0]\n";
			print LOGFILE "Warning: no self hits found for query $fields[0]\n";
		}
	# no lineage found
		unless (defined $fields[13] && !($fields[13] =~ /no\sdata/) )
		{
			print STDERR "Warning: no lineage found for query $fields[0]\n";
			print LOGFILE "Warning: no lineage found for query $fields[0]\n";
		}
	# besthits with no best hit
		unless (defined $fields[4] && !($fields[4] =~ /no\sdata/) )
		{
			print STDERR "Warning: no best hit found for query $fields[0]\n";
			print LOGFILE "Warning: no best hit found for query $fields[0]\n";
		}
	# besthits with no LPI calculation
		unless (defined $fields[5] && !($fields[5] == 0) )
		{
			print STDERR "Warning: no LPI calculated for query $fields[0]\n";
			print LOGFILE "Warning: no LPI calculated for query $fields[0]\n";
		}
	# presence of self identifier in final results
		my $tax_id = $fields[12];
		my $species_name = $fields[13];
		my $lineage = $fields[14];
		foreach my $key (keys %exclude_terms)
		{
			if (($exclude_terms{$key} eq "name")
			&& ($species_name =~ /$key/))
			{
				print STDERR  "Error: found self-term $key for query $fields[0]\n";
				print LOGFILE "Error: found self-term $key for query $fields[0]\n";
				last;
			}
			elsif (($exclude_terms{$key} eq "number"
			&& ($tax_id == $key))
			)
			{
				print STDERR  "Error: found self-tax-id $key for query $fields[0]\n";
 				print LOGFILE "Error: found self-tax_id $key for query $fields[0]\n";
 				last;
			}
		}
	}	
	close INFILE;
	close OUTPUT;
}

sub histogram_calcs
 {
 	my($min, $max, $num_bins, $values) = @_;
 	
 	my %bins = ();
 	my @sorted_values = sort {$a <=> $b} @$values;  #numeric sort	

# get bin range and sizes, put into hash structure for counting	
	my $range = abs($max - $min);	
	my $binsize = $range/$num_bins;	
	my $next_bin = $min;
	$bins{$next_bin} = 0;
	while ($next_bin < $max)
	{	
		$next_bin += $binsize;
		if ($next_bin > $max)
		{
			$next_bin = $max;
			$bins{$next_bin} = 0;
			last;
		}
		$bins{$next_bin} = 0;
	}
	
# count values in each bin
	my @bin_list = sort {$a <=> $b} keys %bins;	
	my $index = 0;
	foreach my $value (@sorted_values)
	{
		chomp $value;
		next unless ($value =~ /[\d\.e-]+/);	#skip non-numeric values
		my $current_bin = $min;
		while ($value > $bin_list[$index])
		{			
			$index++;		
		}
		$current_bin = $bin_list[$index];
		$bins{$current_bin}++;		
	}
	
	return %bins; 
 }
__END__
