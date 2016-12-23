#!/usr/bin/perl
# generate_dh_self_keywords.pl
# Sheila Podell
# January 25, 2016

# Based on recursive SQL queries to NCBI taxdump names and nodes files,
# tries to lists of tax id's for DarkHorse self-exclusion lists 
# reports at genus, species and strain level taxonomic granularity

# Takes as input a tab-delimited list of genome sequences with the following fields:
	# genome name
	# ncbi_tax_id_number (optional - can be left blank or give value of zero)
# lines starting with "#" are ignored
# Program will return a warning if genome not found in NCBI taxdump database
	# Note: NCBI database sometimes contains multiple redundant taxonomy ids and/or names 
	# associated with the exact same genus, species or strain.

# results go to STDOUT in tab delimited format with the following fields:
	# input_genome_name 
	# primary_tax_id 
	# NCBI_species name 
	# genus_self_terms (comma delimited list) 
	# species_self_terms (comma delimited list) 
	# strain_self_terms (comma delimited list) 

use warnings;
use strict;
use Getopt::Long;
use DBI;

# get command line arguments
	my ($db, $datafile, $sth, $sql, $config_filename, $printline, $debug);
	my $USAGE;
	my $message = qq(
	Usage: $0  -i source_datafile -c config filename
	source_datafile format (tab-delimited): 
	species_name\tncbi_tax_id)."\n\n";
	GetOptions( 
			    "i=s" => \$datafile,
			    "c=s" => \$config_filename,
			    "d=s" => \$debug
				);	     
				
# check arguments
	if($USAGE || !$datafile)
	{
		print STDERR $message;
		exit(0);
	} 	
	unless (-s $datafile)
	{
		print STDERR "  Couldn't find source datafile $datafile \n";
		print LOGFILE "  Couldn't find source datafile $datafile \n";
		exit(0);
	}	
	unless (-s $config_filename)
	{
		print STDERR "  Couldn't find config file $config_filename \n";
		print LOGFILE "  Couldn't find config file $config_filename \n";
		exit(0);
	}	

# get SQL database parameters from config file
	my %config_params = &read_config_file($config_filename);
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
	my @row = ();
	my $element;
	
# open a logfile, put STDERR comments there too.
	my $logfilename = "$$"."_logfile";
	my $date = `date`;
	chomp $date;
	open (LOGFILE, ">$logfilename") or die "can't open log file $logfilename\n $!\n";
	print LOGFILE qq($date
$0
input file = $datafile
config file = $config_filename 
);
	
# get list of input names, create genome objects
	my @genome_list = ();		# original user entry order
	my %genome_objects = ();	# key = user-entered genome name, value = object ref
	my %uniq_primary_ids =();
	my $current = "";
    open (INPUT, $datafile)  or die "can't input file $datafile\n$!\n";
    while (<INPUT>)
    {
   		chomp;
   		next if $_ =~ /^\s*$/;	# ignore blank lines
   		next if $_ =~ /^#/;	    # ignore comments
   		next if $_ =~ /^genome/;	# ignore header line
   		my @tmp = split "\t", $_; 		
   		unless (scalar @tmp >0)
   		{
   			print STDERR "Input data error, line $.:\n$_\n";
   			print LOGFILE "Input data error, line $.:\n$_\n";
   			print STDERR "Expected format: species name<tab>ncbi_tax_id \n";
   			next;
   		}  		
   	# remove extra spaces in or after genome names
   		my $genome_name = $tmp[0];
		$genome_name =~ tr/ / /s;
		if ($genome_name =~ /(.+)\s+$/)
		{
			$genome_name = $1;
		}
   		if (exists $genome_objects{$genome_name})
		{
			print STDERR "  WARNING: input contains duplicate genome name: $_\n";
			print LOGFILE "  WARNING: input contains duplicate genome name: $_\n";
			next;
		}
	
	# exclude genomes missing tax_id, genus, and species name - too difficult to get from database 
		unless (defined $tmp[1] && $tmp[1] >0)
		{
			if ($genome_name =~ /unclassified/ || $genome_name =~ /unidentified/)
			{
				print STDERR "  ERROR: skipping genome \"$genome_name\" - can't be classified.\n";
				print LOGFILE "  ERROR: skipping genome \"$genome_name\" - can't be classified.\n";
				next;
			}
		}		

	# create object to hold genome attributes	
		my $primary_tax_id =  $tmp[1];
		if (defined $primary_tax_id && $primary_tax_id > 0)
		{
			$uniq_primary_ids{$primary_tax_id}++;
			if ($uniq_primary_ids{$primary_tax_id} > 1)
			{
				print STDERR "  ERROR: non-unique primary tax id $primary_tax_id, $genome_name";
				print LOGFILE "  ERROR: non-unique primary tax id $primary_tax_id, $genome_name";
			}
		}
		my @args = ($genome_name, $primary_tax_id);
		my $current = &new_genome("genome", \@args);
		if (defined $current) 
		{
			$genome_objects{$genome_name} = $current;
			push @genome_list, $genome_name;
		}	
	}
  close INPUT;
  
# Count number of genomes, give feedback
	my $num_genomes = scalar keys %genome_objects;
	my $num_analyzed = 0; 
	unless ($num_genomes >0)
	{
		print STDERR "ERROR: no valid data to process\n";
		print LOGFILE "ERROR: no valid data to process\n";
		exit(0);
	}	  
	print STDERR "\nFound $num_genomes genomes to analyze . . .\n";
	print LOGFILE "\nFound $num_genomes genomes to analyze . . .\n";
	
# start output
	my @col_headers = ("genome_name", "primary_tax_id", "ncbi_name", 
			"genus_list", "species_list", "strain_list");
	my $headers = join "\t", @col_headers;
	print "$headers\n";
 
# get tax_id and name info from SQL databases 
foreach my $genome_name (@genome_list)
{
	my $genome_obj = $genome_objects{$genome_name};
	my $primary_tax_id = $genome_obj->{primary_tax_id};
	my $current_id;
	my $current_name = $genome_name;
	my $rank = "";
	my $ncbi_name = "ncbi_name";
	my $genus_name = "";
	my $counter = 0;  # variable to prevent run-away recursion
		
# if no primary tax id given, try to find it using genome name
	unless (defined $primary_tax_id && $primary_tax_id > 0)
	{
		if ($debug) {print STDERR "\nSearching for tax id, $genome_name\n";}
		$primary_tax_id = &get_primary_id($genome_name);
		unless (defined $primary_tax_id && $primary_tax_id > 0) 
		{
			print STDERR "  WARNING: no primary tax_id found for $genome_name";
			print LOGFILE "  WARNING: no primary tax_id found for $genome_name";
			$primary_tax_id = 0;
			next;
		}
		$genome_obj->{primary_tax_id} = $primary_tax_id;		
	}

	# get query-specific tax ids and names by traversing linked lists in ncbi nodes table 
	if (defined $primary_tax_id && $primary_tax_id > 0)
	{							
		$counter = 0;
		$current_id = $primary_tax_id;
		
		while ($counter < 4)
		{				
			$counter++;
			
			my $result_ref = &get_tax_id_rank($current_id);
			my @result = @$result_ref;
			
			my $new_rank =  $result[0];
			my $new_name =  $result[1];
			my $new_name_class =  $result[2];
			
			# remove weird characters from name?
			
					
			if (defined $new_rank && length $new_name > 3
				&& defined $new_rank && length $new_rank >3)
			{
				$rank = $new_rank;
				$current_name = $new_name;
			}
			else 			
			{
				if ($debug) {print STDERR "  ERROR: undefined name/rank,input #$num_analyzed, $genome_name\n";}
				print LOGFILE "  ERROR: undefined name/rank, input #$num_analyzed, $genome_name\n";
				last;
			}
							
			if ($rank eq "no rank" || $rank eq "subspecies"
			&& scalar @{$genome_obj->{subspecies}} < 2)
			{
				push @{$genome_obj->{subspecies_list}}, $current_id;
				$genome_obj->{strain_name} = $current_name;
				$genome_obj->{strain_id} = $current_id;			
			}
			if ($rank eq "species")
			{
				push @{$genome_obj->{species_list}}, $current_id; 
				$genome_obj->{species_name} = $current_name;
				$genome_obj->{species_id} = $current_id;	
			}
			if ($rank eq "genus")
			{
				push @{$genome_obj->{genus_list}}, $current_id;
				# get $ncbi_name from the tax_id, put into object.
					$genome_obj->{genus_name} = $current_name;
					$genome_obj->{genus_id} = $current_id;
					if ($current_name =~ /^Candidatus\s(.*)/)
					{
						$current_name = $1;
					}
					push @{$genome_obj->{genus_list}}, $current_name;;
			}
			
			if ($counter == 1)
			{
				if ($debug) {print STDERR "Analyzing primary id $current_id for $genome_name\n";}
				print LOGFILE "Analyzing primary id $current_id for $genome_name\n";
				#$genome_obj->{ncbi_name} = $current_name;
			}
					
			last if ($genome_name =~ /uncultured/ || $genome_name =~ /unidentified/
				|| $genome_name =~ /marine/);
			my $next_id = &get_parent_tax_id($current_id);
			$current_id = $next_id;		
			last if ($rank eq "genus");
			
		}
		$num_analyzed++;
		if ($num_analyzed % 10 == 0)
		{
				print STDERR "\t$num_analyzed genomes analyzed so far.\n";
				print LOGFILE "## $num_analyzed genomes analyzed so far.\n";
		}
		&write_output($genome_name);
		
	}
	else 
	{
		print STDERR "  WARNING: no primary tax id found for $genome_name\n";
		print LOGFILE "  WARNING: no primary tax id found for $genome_name\n";
		next;
	}					
}

print STDERR "  Finished analyzing $num_analyzed genomes\n";
$date = `date`;
print LOGFILE "  Finished analyzing $num_analyzed genomes\n$date\n";

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

sub new_genome
{
 my ($className, $param) = @_;
  my $self = {};
  bless $self, $className;
  my @properties = @$param;  

# Constructor attributes 
	$self->{genome_name} = "$properties[0]";
	$self->{primary_tax_id} = "$properties[1]" || "0";
	
	$self->{ncbi_scientific_name} = "";
	my @genus_list = ();
	my @species_list = ();
	my @strain_list = ();
	$self->{genus_list} = \@genus_list;
	$self->{species_list} = \@species_list;
	$self->{strain_list} = \@strain_list;
	
   return($self)
}

sub get_primary_id
{
	my ($name) = @_;
	my @terms = split " ", $name;
	my $output;
	print STDERR "Finding primary id for $name\n";
	
	# hacks for bad genus names
			if ($terms[0] eq "marine" && $terms[1] eq "gamma" && $terms[2] eq "proteobacteria")
			{
				$terms[0] = "marine gamma proteobacteria";
			
			}						
			if ($terms[0] eq "Candidatus")
			{
				my $compound_name = "$terms[0]"." $terms[1]";
				$terms[1] = $compound_name;
				shift @terms;
			}
			if ($terms[0] eq "uncultured" && $terms[1] eq "Termite" 
				&& $terms[2] eq "group" && $terms[3] eq "group")
			{
				my $compound_name = "";
				$terms[1] = $compound_name;
				shift @terms;		
			}	
			
			if ($terms[0] eq "unidentified" && $terms[1] eq "eubacterium")
			{
				my $compound_name = "$terms[0]"." $terms[1]";
				$terms[1] = $compound_name;
				shift @terms		
			} 
			if ($terms[0] eq "uncultured" && $terms[1] eq "bacterium")
			{
				my $compound_name = join " ", @terms;
				@terms = ($compound_name);		
			} 		 
	
	# try finding names matching the one given by progressively trimming off the final term
	for (my $i=0; $i< scalar @terms; $i++)
	{
		my @slice = @terms[0..$i];			
		my $query = join " ", @slice;
		
		if ($i == scalar @terms -1)
		{
			
			$sql = qq(select distinct tax_id from names where name_txt = '$query');
		}
		else
		{
			$sql = qq(select distinct tax_id from names where name_txt = '$query');
		}
		
		$sth = $dbh->prepare($sql)
				or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
		$sth->execute()
				or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");
			
			@row = ();
			my $output;
			my $id_count = 0;
			while (@row = $sth->fetchrow_array)
			{
		 		next if $row[0] =~ /tax_id/;	#don't include SQL header
				$output .= "$row[0]";
				$id_count++;
			}
		if (defined $output && length $output >1)
		{
			# remove terminal or extra space characters
			$output =~ tr/ / /s;
			if ($output =~ /(.+)\s+$/)
			{
				$output = $1;
			}
			return $output;
		}			
	}
}

sub get_tax_id_rank
{			
	my ($tax_id) = @_;
	my $current_rank;
	my $current_name;
	$sql = "select rank, name_txt, name_class
	from names, nodes
	where names.tax_id = nodes.tax_id
	and names.tax_id = $tax_id;";
	
	if ($debug) {print STDERR "\n$sql\n";}		
		$sth = $dbh->prepare($sql)
			or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
		$sth->execute()
			or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");

		@row = ();
		while (@row = $sth->fetchrow_array)
		{
			next if $row[0] =~ /rank/;	#don't include SQL header
			next if $row[2] =~ /type\smaterial/;
			$current_rank = "$row[0]";
			$current_name = "$row[1]";
			last if ("$row[2]" eq "scientific name");
		}
		my @out_array = ($current_rank,$current_name);
		my $out_ref = \@out_array;	
		return $out_ref;	 
	}

sub get_parent_tax_id
{			
	my ($input_tax_id) = @_;
	my $parent = "not_found"; 
	$sql = "select parent_tax_id from nodes where tax_id = $input_tax_id;";
	if ($debug) {print STDERR "\n$sql\n";}	
	
	$sth = $dbh->prepare($sql)
		or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
	$sth->execute()
		or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");

	@row = ();
	my $id_count = 0;
	while (@row = $sth->fetchrow_array)
	{
		next if $row[0] =~ /tax_id/;	#don't include SQL header
		$parent = "$row[0]";
	}	
	#if ($debug) {print STDERR "$primary_tax_id $rank $parent_tax_id\n";}
	return $parent;	 
}

sub get_siblings
{		
		my ($parent_id) = @_;
		my @sibs_list = ();
		
		unless (defined $parent_id && $parent_id > 0)
		{
			return undef;
		}
		$sql = "select tax_id  
from nodes 
where parent_tax_id = $parent_id;";
if ($debug) {print STDERR "\n$sql\n";}	

	$sth = $dbh->prepare($sql)
		or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
	$sth->execute()
		or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");

	@row = ();
	my $id_count = 0;
	while (@row = $sth->fetchrow_array)
	{ 
		next if $row[0] =~ /tax_id/;	#don't include SQL header
		push @sibs_list, $row[0];
	}			
		return \@sibs_list;
}

sub write_output
{
	my ($key) = @_;
	my ($genus_list, $species_list, $strain_list);
		my $genome_name = $genome_objects{$key}->{genome_name};
		my $primary_tax_id = $genome_objects{$key}->{primary_tax_id};
		my $ncbi_name = $genome_objects{$key}->{ncbi_name} || "not found";
		
		my @genus_array = @{$genome_objects{$key}->{genus_list}};
		my @species_array = @{$genome_objects{$key}->{species_list}};
		my @strain_array = @{$genome_objects{$key}->{strain_list}};
		
		# make sure none of the lists are empty 
			if (scalar @genus_array ==0)
			{
				push @genus_array, $genome_objects{$key}->{primary_tax_id};
			}
			if (scalar @species_array ==0)
			{
				push @species_array, $genome_objects{$key}->{primary_tax_id};
			}
			if (scalar @strain_array ==0)
			{
				push @strain_array, $genome_objects{$key}->{primary_tax_id};
			}
		
		# get ids for taxonomically related strains (subspecies or no rank with same parent as query species) 
	my $strain_sibs_ref = &get_siblings($genome_objects{$key}->{species_id});
	my @strain_sibs;
	if (defined $strain_sibs_ref)
	{ 
		@strain_sibs= @$strain_sibs_ref;
		push @strain_array, @strain_sibs;
		my $num_sibs = scalar @strain_array;
		#print STDERR "    $num_sibs strain sibs for $genome_name\n";
		print LOGFILE "    $num_sibs strain sibs for $genome_name\n";
	}	
	# get ids for taxonomically related species (same parent genus as query genome)
		my $species_sibs_ref = &get_siblings($genome_objects{$key}->{genus_id});
		my @species_sibs;
		if (defined $species_sibs_ref)
		{ 
			@species_sibs = @$species_sibs_ref;
			push @species_array, @species_sibs;
			my $num_sibs = scalar @species_array;
			print LOGFILE "    $num_sibs species sibs for $genome_name\n";
		}
		print LOGFILE "\n";
	# remove internal replicates from arrays?		
		$genus_list = join ",", @genus_array;
		$species_list = join ",", @species_array;
		$strain_list = join ",", @strain_array;
			
		my @attributes = ($genome_name, $primary_tax_id, $ncbi_name, 
			$genus_list, $species_list, $strain_list); 
		
		my $printstring = join "\t", @attributes;
		print  "$printstring\n";	
	}


__END__

	
############
mysql> select distinct name_class from names;
+---------------------+
| name_class          |
+---------------------+
| synonym             |
| scientific name     |
| in-part             |
| blast name          |
| genbank common name |
| authority           |
| misspelling         |
| type material       |
| equivalent name     |
| includes            |
| genbank synonym     |
| common name         |
| misnomer            |
| genbank acronym     |
| anamorph            |
| genbank anamorph    |
| teleomorph          |
| acronym             |
+---------------------+
18 rows in set (1.87 sec)


mysql> select distinct rank from nodes;
+------------------+
| rank             |
+------------------+
| no rank          |
| superkingdom     |
| genus            |
| species          |
| order            |
| family           |
| subspecies       |
| subfamily        |
| tribe            |
| phylum           |
| class            |
| forma            |
| suborder         |
| superclass       |
| subclass         |
| varietas         |
| kingdom          |
| subphylum        |
| superfamily      |
| infraorder       |
| infraclass       |
| superorder       |
| subgenus         |
| parvorder        |
| superphylum      |
| species group    |
| species subgroup |
| subtribe         |
| subkingdom       |
+------------------+
29 rows in set (0.68 sec)

	