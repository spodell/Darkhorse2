#!/usr/bin/perl
# install_darkhorse2.pl 
# Sheila Podell 
# December 22, 2016

# Tests availability of essential pre-requisites for running DarkHorse
# Generates customized DarkHorse config file based on template and user input
# Automates process of installing MySQL databases required to run DarkHorse
# Generates customized reference protein fasta file to be used in blast (diamond) searches

use warnings; 
use strict;
use Getopt::Long;
use Cwd;
use File::Basename;
use File::Spec;
use DBI;

# global input:
    my ($program_directory, $input_directory);
    my ($db, $sth, $sql, $overwrite, $update, $result);
    my ($table_name, $firstline, $datafile);
    my $config_template = "templates/config_template";
    my $config_filename = "$config_template";
    my $user_message = "";
    my $cmd = "";
    my $debug = 0;
    my $num_cpu = 1;
    
    my $USAGE;
	my $message = qq(
  Usage: $0  -p program directory path 
  Optional parameters: 
    -c configuration file
    -i input data directory
    -o allow database overwrites
    -u allow database updates)."\n\n";
    
	GetOptions( 
			    "p=s" => \$program_directory,
			    "i=s" => \$input_directory,
			    "c=s" => \$config_filename,   # pre-existing or manually edited file
			    'o'   => \$overwrite,		
			    'u'   => \$update,
			    "n=i" => \$num_cpu,
			    "d=i" => \$debug,
				);
	 &check_args();

# Get default parameters from template, then generate customized config file based on user input    
    my %config_params = &read_config_file("$program_directory/$config_template");
    my $new_filename = &new_filename();
    if ($config_filename ne "$config_template" && -s "$program_directory/$config_filename")
    {
        print "  Using config file $config_filename\n";
    }     
    else
    {
        print "  Creating customized configuration file: $new_filename\n";
        $config_params{program_directory} = $program_directory;
        &write_config_file(\%config_params, $new_filename);
        $config_filename = $new_filename;
    }
  # revise parameters interactively if necessary
     %config_params= &read_config_file($config_filename);
     &user_input(\%config_params);
     
# set performance monitoring variables
	my ($date, $subroutine_multiply_factor, $est_time_mins, $est_time_hours, $text_with_commas); 
    my $sql_time_factor = ($config_params{max_lines_per_packet})/2000;	
	my $max_lines_per_hash = $config_params{max_lines_per_packet} * 1000;  # base level 2 million

# determine absolute paths for logfile
	my $abs_prog_directory = File::Spec->rel2abs( $program_directory);
	my $abs_config_filename = File::Spec->rel2abs("$program_directory/$config_filename");
		
# create logfile to keep track of progress
	my $logfile = "darkhorse_install_"."$$".".log";
	open (LOGFILE, ">$logfile") or die "couldn't open logfile $logfile for writing, $!\n";	
	print LOGFILE qq(###################################################
  Darkhorse version 2.0 (beta) installation logfile
###################################################

Installation settings: 
	program directory = $abs_prog_directory
	config filename = $abs_config_filename
	reference fasta file = $config_params{genbank_nr_fasta_path}
	accession tax id file = $config_params{protid_taxid_dmp_path}
	MySQL database name = $config_params{db_name}
	max_lines_per_packet = $config_params{max_lines_per_packet}
			    
);
  
# Install downloaded NCBI taxonomy names table in local MySQL database
    my $install_script_dir = "$program_directory/bin/installation_scripts";
    my $database = $config_params{db_name};
    my $names = $config_params{names_dmp_path};
    my $nodes = $config_params{nodes_dmp_path};
    my $accession_taxid_index = $config_params{protid_taxid_dmp_path};
    my $nr_fastapath = $config_params{genbank_nr_fasta_path};
    $date = `date`;
    chomp $date;
    $user_message = "$date"."  Loading taxonomy names table into MySQL database $database";
    &send_message("$user_message");	
   $cmd= "$install_script_dir/load_taxonomy_names.pl -d $database -c $config_filename -i $names";
   `$cmd`;
   my $num_name_rows = &check_table_size(\%config_params,"names");
   
# Install downloaded NCBI taxonomy nodes table in local MySQL database
   $date = `date`;
   chomp $date;
   $user_message =  "$date"."  Loading taxonomy nodes table into MySQL database $database";
    &send_message("$user_message");
   $cmd = "$install_script_dir/load_taxonomy_nodes.pl -d $database -c $config_filename -i $nodes";
   `$cmd`;    
	my $num_node_rows = &check_table_size(\%config_params,"nodes");
	
# get list of unique tax id numbers from prot.accession2taxid
	my($accession_taxid_index_shortname, $dirs, $suffix) = fileparse($accession_taxid_index);
	$date = `date`;
	chomp $date;
	$user_message = "$date"."  Collecting unique tax ids from $accession_taxid_index_shortname";
	&send_message($user_message);

# estimate time based on size of input file
	my $accession_taxid_index_filesize = &check_filesize($accession_taxid_index);
		
# sort accession index on taxid, create file for use in multiple subsequent steps
	$cmd = "LC_ALL=C sort -n -k3 $accession_taxid_index | cut -f2,3 > $input_directory/sorted_accession_taxid_index";
	`$cmd`;
	$cmd = "cut -f2 $input_directory/sorted_accession_taxid_index | uniq | LC_ALL=C grep -v taxid > $input_directory/catalog_tax_ids";
	`$cmd`;
	my $num_prelim_tax_ids = 0;
	my $count = `wc -l $input_directory/catalog_tax_ids`;
	chomp $count; 
		
	if ($count =~ /(\d+)/)
	{
		$num_prelim_tax_ids = $1;
	}
		
	unless ($num_prelim_tax_ids >0)
	{
		$user_message = "ERROR"."  Found $num_prelim_tax_ids unique taxonomy ids";
		&send_message($user_message);
		exit(0);
	}
	$date = `date`;
	chomp $date;
	$text_with_commas = &commify($num_prelim_tax_ids);
	$user_message = "$date"."  Found $text_with_commas unique taxonomy ids.";
	&send_message($user_message);
	
# pre-calculate lineages by recursively traversing nodes from NCBI taxonomy database
# note: recursive_taxonomy_taxid.pl skips tax_ids 0,1,2 (all, root, Bacteria) 
# skips obsolete tax_ids with warning (e.g. prot.accession2taxid not synchronized with taxonomy names table)
  
	$date = `date`;
	chomp $date;
    $user_message = "$date"."  Pre-calculating lineages.";
    $est_time_mins = int($num_prelim_tax_ids * 0.0003)/$sql_time_factor;
    if ($est_time_mins > 0)
    {
		&estimate_time($est_time_mins, $user_message);
	}
	 	
  $cmd = "$install_script_dir/recursive_taxonomy_taxid.pl -d $database -c $config_filename -i $input_directory/catalog_tax_ids > $input_directory/catalog_tax_id_lineages";	
	`$cmd`;

   unless (-s "$input_directory/catalog_tax_id_lineages")
 	{
 		$user_message = "ERROR: Creation of file $input_directory/catalog_tax_id_lineages unsuccessful\n";
 		&send_message($user_message);
 		exit(0);
 	}
 	
 	unless ($debug > 0)
 	{
 		&unlink_file("$input_directory/catalog_tax_ids");
 	}
 	
 # get species names for tax_ids - needed for precalc_lineages table loading
	$est_time_mins = int($num_prelim_tax_ids * 0.00001)/$sql_time_factor;
	$est_time_mins = sprintf "%01X", $est_time_mins;
	
	$date = `date`;
	chomp $date;
    $user_message = "$date"."  Collecting species names associated with tax_ids.";
    if ($est_time_mins > 1)
    {
    	$user_message .= "\n   NOTE: this step may take $est_time_mins or more minutes.";
    }   
    &send_message($user_message);
    
# conserve memory by writing species names directly to external file      
    &get_species_names("$input_directory/catalog_tax_id_lineages","$input_directory/raw_lineages.tab");
    
      unless (-s "$input_directory/raw_lineages.tab")
 	{
 		$user_message = "ERROR: Creation of file $input_directory/raw_lineages.tab unsuccessful\n";
 		&send_message($user_message);
 		exit(0);
 	}
 	
	 unless ($debug > 0)
	 {
 		&unlink_file("$input_directory/catalog_tax_id_lineages");
	 }
   	           
# load precalc_lineages table
	$date = `date`;
	chomp $date;
    $user_message = "$date"."  Loading precalc_lineages table into MySQL database $database";
    &send_message($user_message);
    $cmd = "$install_script_dir/load_precalc_lineages.pl -d $database -c $config_filename -i $input_directory/raw_lineages.tab";     
	`$cmd`;
	my $num_precalc_lineages_rows = &check_table_size(\%config_params,"precalc_lineages");
	my $lineage_diff = $num_prelim_tax_ids - $num_precalc_lineages_rows;
	if ($lineage_diff >0)
	{
		#$user_message = "\n  NOTE: In $lineage_diff cases, more than one tax_id was associated with the same lineage (this is OK).";
		#&send_message($user_message);
	}
	
	unless ($debug >0)
	{		
		&unlink_file("$input_directory/raw_lineages.tab");
	}
	
# create list of valid, informative tax ids (exclude uninformative lineages, )
	
	my $informative_taxid_filename = "$input_directory/informative_taxids";
	$date = `date`;
	chomp ($date); 
	$user_message =  "$date"."  Identifying informative tax_ids";
	&send_message($user_message);			
	&select_informative_taxids($informative_taxid_filename);
			 
# load prot_tax_id table, including only sequences informative tax_ids
	$date = `date`;
	chomp ($date); 
	$user_message =  "$date"."  Loading table protid_taxid into MySQL database $database.";
	&send_message($user_message);
	 
    $est_time_mins = int($accession_taxid_index_filesize * 0.0005)/$sql_time_factor;
    if ($est_time_mins > 0)
    {
		&estimate_time($est_time_mins, $user_message);
	}

	 $cmd = "$install_script_dir/load_protid_taxid_table.pl -c $config_filename -i $input_directory/sorted_accession_taxid_index -v $informative_taxid_filename";
	 `$cmd`; 
	$date = `date`;
	chomp ($date);
	$user_message = "$date"."  protid_taxid table loading complete.";
	my $num_protid_taxid_rows = &check_table_size(\%config_params,"protid_taxid");
	unless ($debug > 0)
	{
		&unlink_file("$input_directory/sorted_accession_taxid_index");
		&unlink_file("$input_directory/informative_taxids"); #$informative_taxid_filename"		
	}	
		
# Estimate time required to process Genbank nr fasta file (annotation, sequence length)
# counting the all the fastas in the current nr file with unix (96,320,599) took 22 minutes			
	my $fasta_file = "$config_params{genbank_nr_fasta_path}"; # or subfile?
	my $fasta_size = (stat $fasta_file)[7];
	$fasta_size = int($fasta_size/1000000); # MB
	
	$date = `date`;
	chomp $date;
    $user_message = "$date"."  Counting reference fasta sequences.";
	
	$est_time_mins = int($fasta_size * 0.0000005); # expect ~1 min per 2 GB
    if ($est_time_mins > 0)
    {
		&estimate_time($est_time_mins, $user_message);
	}
		
	$cmd = "LC_ALL=C grep '>' $fasta_file |wc -l" ; 
		my $num_ref_sequences = `$cmd`;
		chomp $num_ref_sequences;
		$num_ref_sequences = int($num_ref_sequences);
		my $num_seqs_per_subfile = int($num_ref_sequences/$num_cpu) + 1;	
		unless ($num_ref_sequences > 0)
		{
				$user_message = "ERROR: no valid reference sequences found.";
				&send_message($user_message);
				exit (0);
		}
		$date = `date`;
		chomp $date;
		$text_with_commas = &commify($num_ref_sequences);
		$user_message = "$date"."  Extracting annotations for $text_with_commas reference sequences.";
		
	$est_time_mins = int($num_ref_sequences * 0.0000007)/$num_cpu;
	
	if ($est_time_mins > 0)
	{
		&estimate_time($est_time_mins, $user_message);
	}
	else
	{
		&send_message($user_message);
	}	
 	  
 # get the fasta info into a file suitable for database loading
 	 $cmd = "$install_script_dir/get_nr_seq_info.pl $nr_fastapath > $input_directory/annotations.tab";	 
 	`$cmd`;
 	unless (-s "$input_directory/annotations.tab")
  	{
  		print STDERR "ERROR: Creation of file $input_directory/annotations.tab unsuccessful\n";
  		exit(0);
  	}	
  	my $num_annot =`wc -l $input_directory/annotations.tab`;
  	chomp $num_annot;
  	if($num_annot =~ /(\d+)/)
  	{
  		$num_annot = $1;
  	}
  	$text_with_commas = &commify($num_annot);
  	$date = `date`;
  	chomp $date;
  	$user_message = "$date"."  Found $text_with_commas annotations.";
  	&send_message($user_message);	

# load annotations table
	$date = `date`;
  	chomp $date;	
	$user_message = "$date"."  Loading annotations table into MySQL database $database.";
	&send_message($user_message);
	$cmd = "$install_script_dir/load_annotation_table.pl -c $config_filename -i $input_directory/annotations.tab";
	`$cmd`;

	$date = `date`;
	chomp ($date);
	$user_message = "$date"."  Annotation table loading complete.";
	my $num_annotation_rows = &check_table_size(\%config_params,"annotation");
	
	unless ($debug > 0)
	{
		&unlink_file("$input_directory/annotations.tab");
	}
			
	$date = `date`;
	chomp ($date);
	$user_message = "$date"."  Combining taxonomy info with sequence data.";
		
	$est_time_mins = int($num_ref_sequences * 0.0000008)/$num_cpu;
	
	if ($est_time_mins > 0)
	{
		&estimate_time($est_time_mins, $user_message);
	}
	else
	{
		&send_message($user_message);
	}

# merge annotation table parse_prot_accession2taxid table
	$cmd = "$install_script_dir/merge_annot_protid_taxid_tables.pl -c $config_filename > $input_directory/nr_all_ids.tab";	
	`$cmd`;
	
	$date = `date`;
	chomp ($date);	
	if (-s "$input_directory/nr_all_ids.tab")
	{
 		$user_message = "$date"."  Finished combining taxonomy info with sequence data";
 		&send_message($user_message);
 	}	
	else 
	{
 		$user_message = "ERROR: Creation of file $input_directory/nr_all_ids.tab unsuccessful\n";
 		&send_message($user_message);
 		exit(0);
 	}
 	
 	&drop_table("protid_taxid");
 
# load db_seq_ids table	
	$date = `date`;
	chomp ($date); 
	 $user_message =  "$date"."  Loading table db_seq_ids into MySQL database $database.";
	 &send_message($user_message);
	 
	 $cmd = "$install_script_dir/load_db_seq_ids_table.pl -d $database -c $config_filename -i $input_directory/nr_all_ids.tab";
	 `$cmd`; 
	$date = `date`;
	chomp ($date);
	$user_message = "$date"."  db_seq_ids table loading complete.";
	my $num_load_db_seq_ids_rows = &check_table_size(\%config_params,"db_seq_ids");	
	unless ($debug >0 )
	{
		&unlink_file("$input_directory/nr_all_ids.tab");
	}
 
# create lineage_index table to speed searches 
    $date = `date`;
	chomp ($date);
	$user_message = "$date"."  Creating master lineage index";
	&send_message($user_message);
    
    `$install_script_dir/collect_lineage_index_data.pl -d $database -c $config_filename -a > $input_directory/all_lineages.tab`;
    
    if (-s "$input_directory/all_lineages.tab")
 	{
 		#$user_message = "$date"."  master lineage index complete.";
 		#&send_message($user_message);
 	}  
    else
 	{
 		$user_message = "ERROR: Creation of file $input_directory/all_lineages.tab unsuccessful.";
 		&send_message($user_message);
 		exit(0);
 	}
	
# load lineage index table
	$date = `date`;
	chomp ($date); 
	 $user_message =  "$date"."  Loading table lineage_index into MySQL database $database.";
	 &send_message($user_message);	  
    `$install_script_dir/load_lineage_index_table.pl  -d $database -c $config_filename -i $input_directory/all_lineages.tab`;

	$date = `date`;
	chomp ($date);
	$user_message = "$date"."  db_seq_ids table loading complete.";
	my $num_load_lineage_index_rows = &check_table_size(\%config_params,"lineage_index");
	
	unless ($debug >0)
	{	
		&drop_table("annotation");
	}		

# user feedback
	$date = `date`;
	chomp ($date); 
	 $user_message =  "$date"."  Finished loading MySQL tables.";
	 &send_message($user_message);
	 
 # create reference fasta file for blast
 	$date = `date`;
	chomp ($date); 
	 $user_message =  "$date"."  Creating informative_ref_seqs.fasta file.";
	 &send_message($user_message);
 	
 	#note: this may be fewer sequences than used to create db_seq_ids, because
 	# some may have had obsolete tax_ids, preventing automated lineage calculation
 	
 	$cmd = "cut -f1 $input_directory/all_lineages.tab > $input_directory/selected_seq_ids";
 	`$cmd`;
 	my $revised_nr_fasta = "$database"."_informative_ref_seqs.fasta";
 	
 	$cmd = "$install_script_dir/getseq_multiple_nr.pl $input_directory/selected_seq_ids $nr_fastapath > $revised_nr_fasta";
 	`$cmd`;
 	
	if (-s "$revised_nr_fasta")
	{
 		$user_message = "$date"."  Finished selecting informative reference sequences into fasta file";
 		&send_message($user_message);
 	}	
	else 
	{
 		$user_message = "ERROR: Creation of informative_ref_seqs.fasta file unsuccessful\n";
 		&send_message($user_message);
 		exit(0);
 	}
	
# Final clean up 
	unless ($debug > 0)
	{
		&unlink_file("$input_directory/all_lineages.tab");
		&unlink_file("$input_directory/selected_seq_ids");		
		&send_message("   deleting temporary directory $input_directory.");
		rmdir "$input_directory";
	}
	
	$date = `date`;
	chomp $date;
    $user_message = "$date"."  Finished installing DarkHorse2. Thanks for your patience.";
    &send_message("$user_message");				
	
close LOGFILE;

#########################################################
# SUBROUTINES
#########################################################
sub send_message
{
	my ($msg, $error) = @_;
	
	print LOGFILE "$msg\n";
	print STDERR "$msg\n";	
}
	
sub check_args
{
 	if($USAGE || !$program_directory) 
	{
		print STDERR $message;
		exit(0);
	} 
	unless (-s "$program_directory")
	{
		print STDERR "  Couldn't find program directory $program_directory \n";
		exit(0);
	}
	unless (-s "$program_directory/$config_template")
	{
		print STDERR "  Couldn't find config template $program_directory/$config_template \n";
		exit(0);		
	}
	# check for required program files
        my @program_list = (
                            "load_taxonomy_names.pl",
                            "load_taxonomy_nodes.pl",
                            "recursive_taxonomy_taxid.pl",
                            "load_precalc_lineages.pl",
                            "load_protid_taxid_table.pl",
                            "get_nr_seq_info.pl",
                            "load_annotation_table.pl",
                            "merge_annot_protid_taxid_tables.pl",
                            "load_db_seq_ids_table.pl",                           
                            "collect_lineage_index_data.pl",
                            "load_lineage_index_table.pl",
                            "getseq_multiple_nr.pl"                           
                            );
         
         foreach (@program_list)
         {
           my $program =  "$program_directory/bin/installation_scripts/$_";
            unless (-s $program && -x $program)
            {
                print STDERR "Couldn't find program $_\n";
                exit (0);
            }         
         }
	# setup data input directory
        if (defined $input_directory && -s "$input_directory")
        {
            return;
        }
            
        $input_directory = "$program_directory/input_data";
        print STDERR "  Using data input directory $program_directory/input_data\n";	
        
        unless (-s 	$input_directory)
        {
            mkdir "$program_directory/input_data", 0755 or die "Cannot make input_data directory $program_directory/input_data\n: $!"; 		
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

sub write_config_file
{  
	my ($hashref, $newfile) = @_;
	my %parameters = %$hashref;
	
	my $date = `date`;
	chomp $date;
	    
	open (CONFIG,  ">$program_directory/$newfile") or die "couldn't create new configuration file\n$!\n\n";			
	print CONFIG "# Configuration file for DarkHorse $newfile, generated by dh_installation.pl\n";
	print CONFIG "# $date\n\n";
				
	foreach (keys %parameters)
	{
	    print CONFIG "[$_]=$parameters{$_}\n";
	
	}		
    print CONFIG "\n";		 
	close CONFIG;
}

sub new_filename
{
    my $newname = "";
    # don't overwrite old file
        if (-s "$program_directory/darkhorse.cfg")
        {	    
            `cp $program_directory/darkhorse.cfg $program_directory/darkhorse_old.cfg`;	
        }
	
	# put a unique, sequential number on the new filename
        my $basename = "darkhorse";
        my $counter = 1;
        do{
            $newname = "$basename$counter".".cfg";
            $counter++;
        }while (-s "$program_directory/$newname" );
        
	return $newname;
}	
	
sub user_input
{    
    my $changes = 0;
   
    # confirm location of non-executable file paths    
        my @file_list = ( "genbank_nr_fasta_path",
                          "names_dmp_path",
                          "nodes_dmp_path",
                          "protid_taxid_dmp_path",
                         );
        print "\n  Please enter corrected location for the following paths, if required\n";
        foreach (@file_list)
        {
            my $path = $config_params{$_};
            $path = " " unless (defined $path);
            
            print "  $_ ($path): \n";
            my $tmp = <STDIN>;
            chomp $tmp;
            if (length $tmp>1)
            {
                $config_params{$_} = $tmp;
                $changes++;
            }
            unless (-s "$config_params{$_}")
            {
                print  "  Can't find $_ in specified location:\n  $config_params{$_}\n";
                print  "  DarkHorse configuration cannot be completed without this file.\n";
                unlink "$program_directory/$new_filename";
                exit (0);      
            }        
        }

    # confirm mysql connection host, username, password, and database
        my @db_params = ("db_name","db_host","db_user","db_user_password");
        print "  Please enter corrected connection information for MySQL database, if necessary.\n";        
        foreach (@db_params)
        {           
            print "  $_ ($config_params{$_}):\n";
            my $tmp = <STDIN>;
            chomp $tmp;
            if (length $tmp>1)
            {
                $config_params{$_} = $tmp;
                $changes++;
            }                   
        }
        my $result = &check_db_connection(\%config_params);
        unless ($result ==1)
        {
            print "\n  Unable to access database using the following parameters:\n";
             foreach (@db_params)
            {           
                print "    $_ ($config_params{$_})\n";    
            }
            print "  If database $config_params{db_name} does not already exist, you must create it manually.\n";
            exit(0);
        }
        print "Database connection OK\n";       
        
    # check format for max_lines_per_packet, minimum number of lineage terms, min_align_coverage
        unless ($config_params{max_lines_per_packet} =~ /^\d+$/)
        {
            print STDERR "Illegal value for max_lines_per_packet:$config_params{max_lines_per_packet}\n";
            exit(0);	
        }
        unless ($config_params{min_lineage_terms} =~ /^\d+$/)
        {
            print STDERR "Illegal value for min_lineage_terms:($config_params{min_lineage_terms}\n";
            exit(0);	
        }
        unless ($config_params{min_align_coverage} >=0 && $config_params{min_align_coverage} <= 1)
        {
            print STDERR "Illegal value for min_align_coverage:($config_params{min_align_coverage}\n";
            exit(0);	
        }
        
        if ($changes>0)
        {
            &write_config_file(\%config_params, $new_filename);
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

sub missing_file_msg
{
       print  "  This program is required for DarkHorse configuration.\n";
                print  qq(  It can be downloaded from ftp://ftp.ncbi.nih.gov/blast/)."\n";
                unlink "$program_directory/$new_filename";
                exit (0);
}

sub estimate_time
{
	my ($time_mins, $user_message) = @_;
	$est_time_mins = int($time_mins);
	$est_time_hours = $time_mins/60;
	$est_time_hours = sprintf "%01X", $est_time_hours;
		
	if ($est_time_hours > 0.99)
	{
		$user_message .= "\n  NOTE: this step may take $est_time_hours or more hours."; 
	}	
	else
	{
		$user_message .= "\n  NOTE: this step may take $est_time_mins or more minutes."; 
	}
	&send_message($user_message);
}
	
sub check_table_size
{
	my ($hashref, $tablename) = @_;
# set up db connection;
	my %dh_config = %$hashref;	
	my $db_name = $dh_config{db_name};
	my $db_program = $dh_config{db_program};
	my $db_host = $dh_config{db_host};
	my $db_user = $dh_config{db_user};
	my $db_user_password = $dh_config{db_user_password};    
    my $db_path = "DBI:$db_program:$db_name:$db_host:";
    my $dbh = DBI->connect($db_path, $db_user, $db_user_password);

	my $sql = "select count(*) from $tablename";
	my $sth = $dbh->prepare($sql)
    	or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
	$sth->execute()
    	or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");
	
	my @row = ();
	my $element = ();
	my $output = 0;
	while (@row = $sth->fetchrow_array)
	{
		$output = $row[0]; 
	}
	
	$dbh->disconnect();
	
	$date = `date`;
	chomp $date;
	
# user output
	if ($output < 1 || $output < 1 )
	{
		$user_message = "$date"."  ERROR loading MySQL table $tablename: contains $output rows";
		&send_message("$user_message");
		exit (0);
	}
	else
	{
		my $text_with_commas = &commify($output);
		$user_message = "$date"."  MySQL $tablename table loading complete: contains $text_with_commas rows";
		&send_message("$user_message");	
	}
	return $output;
}

sub commify {
	# per http://docstore.mik.ua/orelly/perl/cookbook/ch02_18.htm
    my ($original_value) = @_;
    my $text = $original_value;
    $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text;
}

sub select_informative_taxids
{
# get rid of anything that isn't going to be useful for taxonomic classification
# two criteria: num lineage terms and valid species name
	# species name important to exclude uninformative entries

	my ($outfile) = @_;
	my ($species_name, $ncbi_id);
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

$sql  = qq (
select ncbi_tax_id, species_name
from precalc_lineages
where num_terms > 2;		
);
	$sth = $dbh->prepare($sql)
		or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
	$sth->execute()
	or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");
	
	my @row = ();
	my $element = ();
		
	while (@row = $sth->fetchrow_array)
	{									
		my ($ncbi_tax_id, $species_name) = @row;
		next unless (defined $ncbi_tax_id && $ncbi_tax_id >0);
		my @tmp = split " ", $species_name;	
	# get rid of inconvenient punctuation in species name
		$species_name =~ s/\[//g;
		$species_name =~ s/]//g;
		$species_name =~ s/'//g;
	# make sure it really is a species name
		next unless (scalar @tmp >1);		
		next if $species_name =~ /NULL/;
		next if $species_name =~ /cloning/;
		next if $species_name =~ /clone/;
		next if $species_name =~ /vector/;
		next if $species_name =~ /cosmid/;
		next if $species_name =~ /synthetic/;
		next if $species_name =~ /construct/;
		next if $species_name =~ /contaminant/;
		next if $species_name =~ /unclassified/;
		next if $species_name =~ /unidentified/;
		next if $species_name =~ /unknown/;
		next if $species_name =~ /untyped/;
		next if $species_name =~ /unspecified/;
		print OUTPUT "$ncbi_tax_id\t$species_name\n";
	}
	$dbh->disconnect();
	close OUTPUT;
}

sub sql_select_list
{		
	my ($list, $table_name, $firstline) = @_;
	my @select_list = @$list;
	my $total_num = scalar @select_list;
	my $sql = $firstline;
	
	my @results_tbl = ();
	my @row = ();
	
 # set up db connection;
	my $db_name = $config_params{db_name};
	my $db_program = $config_params{db_program};
	my $db_host = $config_params{db_host};
	my $db_user = $config_params{db_user};
	my $db_user_password = $config_params{db_user_password};    
    my $db_path = "DBI:$db_program:$db_name:$db_host:";
    my $dbh = DBI->connect($db_path, $db_user, $db_user_password);

	my $max_lines_per_packet = $config_params{max_lines_per_packet};
		
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

 sub get_species_names
{
	# only want species names that are present in catalog tax ids
	# process in chunks, to avoid filling memory
		
	my ($datafile, $outfile) = @_;
	$table_name = "names";
	$firstline = 'select name_txt, tax_id from names where name_class = "scientific name" and tax_id in (';
		
	my @id_list =();
	my %tmp_hash = (); #key = tax_id   value = lineage;
	open (INPUT,"$datafile") or die "can't open temp file $datafile for reading\n $!\n";
	open (OUTPUT, ">$outfile") or die "can't open temp file $outfile for writing\n $!\n";

	while (<INPUT>)
	{
		chomp;
		next if ($_ =~ /^\s*$/);
		my @tmp = split "\t", $_;
		my $tax_id = $tmp[0];
			next unless ($tax_id =~ /^\d+$/);
		my $lineage = $tmp[1];
			unless (defined $lineage && length $lineage > 0)
			{
				print STDERR "lineage problem line $. $_\n";
				exit(0);
			}
		push @id_list, $tax_id;
		$tmp_hash{$tax_id} = $lineage;
			
		if (($. % $max_lines_per_hash == 0) || ($. == $num_prelim_tax_ids) || eof)
		{
		# use sql query to get data into array
			my @outlines=&sql_select_list(\@id_list,$table_name, $firstline);
		# append output to new file
			foreach my $next_line (@outlines)
			{
				my ($name_txt, $id) = split "\t", $next_line;
				 # get rid of inconvenient punctuation 
					# $name_txt =~ s/\'//g;
# 					$name_txt =~ s/\\//g;
# 					$name_txt =~ s/\;/,/g;
# 					$name_txt =~ s/\"//g;			
# 					$name_txt =~ s/[//;
# 					$name_txt =~ s/]//;
				my $lineage = $tmp_hash{$id};
				my $out_text = "$name_txt"."\t"."$id"."\t"."$lineage";
				print OUTPUT "$out_text\n";
			}
		# feedback to user 
			if ($. % 100000 == 0)
			{
				#print STDERR "  $. lines processed \n";
			}
		# clear array memory
			@id_list = ();
			$tmp_hash{$tax_id} = ();
		}
	}
	close INPUT;
	close OUTPUT;
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
	print STDERR "     input file size = $filesize_txt\n";
}

sub unlink_file
 {
 	my ($unlink_name) = @_;
 	&send_message("   deleting temporary file $unlink_name.");
 	unlink $unlink_name;
 	
 }
 
sub drop_table
{
	my ($drop_table) = @_;
	&send_message("   dropping temporary MySQL table $drop_table.");

# set up db connection;
	my $db_name = $config_params{db_name};
	my $db_program = $config_params{db_program};
	my $db_host = $config_params{db_host};
	my $db_user = $config_params{db_user};
	my $db_user_password = $config_params{db_user_password};    
    my $db_path = "DBI:$db_program:$db_name:$db_host:";
    my $dbh = DBI->connect($db_path, $db_user, $db_user_password);

	$sql = "drop table $drop_table";
	$sth = $dbh->prepare($sql)
		or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
	$sth->execute()
	or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");
    
    $dbh->disconnect();			
}	


__END__
