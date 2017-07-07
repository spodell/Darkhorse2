#!/usr/bin/env perl
# recursive_taxonomy_v2.pl
# Sheila Podell
# December 16, 2016

# recursively parses Taxonomy database to get lineage for a species
# takes as input, tab-delimited file containing list of tax_ids
# outputs semicolon-delimited lineage for each entry to STDOUT:

	# tax_id lineage

# note: this version skips uninformative input tax_ids 0,1,2 (all, root, Bacteria) 
# warns & skips obsolete tax_ids (e.g. prot.accession2taxid not synchronized with taxonomy names table)
  
use warnings;
use Getopt::Long;
use DBI;

# global variables
	my $recursion_max = 35;
	my $recursion_count = 0;
	
# get command line arguments, defaults 
	my ($db, $datafile, $sth, $sql, $config_filename);
	my $USAGE;
	my $message = qq(
	Usage: $0 -d database_name -i source_datafile -c config filename <-o allow overwrites> )."\n\n";
	GetOptions( "d=s" => \$db,
			    "i=s" => \$datafile,
			    "c=s" => \$config_filename,
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
	
# connect to database
	my $db_path = "DBI:$db_program:$db_name:$db_host:";
	my $dbh = DBI->connect($db_path, $db_user, $db_user_password) or die "Cannot connect: $DBI::errstr"; 	

# set up output order for hierarchical keys - stop after genus
	my @output_list = (
		"superkingdom",
		"kingdom",
		"superphylum",
		"phylum",
		"no rank",
		"subphylum",
		"superclass",
		"class",
		"subclass",
		"infraclass",
		"superorder",
		"order",
		"parvorder",
		"suborder",
		"infraorder",
		"superfamily",
		"family",
		"subfamily",
		"tribe",
		"subtribe",
		"genus",
		"subgenus",
		#"species group",
		#"species subgroup",
		#"species",
		#"subspecies",
		#"forma",
		#"varietas",
		);
	my %output_order =();
	my $index = 0;
	foreach (@output_list)
	{
		$output_order{$_} = $index;
		$index++;
	}
	
# process data	
	open(INPUT,"$datafile")|| die "cannot find input data file $datafile\n";
	while($line=<INPUT>)
	{
        my @rows= ();
        my $num = 0;
        my %query = ();
        my %Tax = (); #key = order retrieved, value = tax_id
        next if($line=~/^#/);
        chomp $line;
        my $input_tax_id = $line;
    # don't process missing or uninformative $input_tax_id (e.g. 0=all, 1=root, 2=Bacteria)
        unless (defined $input_tax_id && length $input_tax_id > 0 && $input_tax_id > 2)
        {
        	#print STDERR "bad input tax id line $.,$input_tax_id.\n";
        	next;
        }        
	
# make sure that the tax id is present in the current darkhorse mysql database
 	next unless (defined $input_tax_id && length $input_tax_id >0);
	$sql = "select tax_id from names where tax_id = $input_tax_id";
	$sth = $dbh->prepare($sql)
		or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
	$sth->execute()
		or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");
	$num = $sth->rows;	
		
# get rest of lineage by recursion, or give up and go on to next line
	if ($num!=0)
	{
		%query = %{$sth->fetchrow_hashref};			
		$tax_id=$query{tax_id};
		%Tax=&nodes($tax_id,%Tax);
	}
# output results in reverse recursion order, removing extraneous stuff		
	unless (defined $query{tax_id})
	{
		print STDERR "  WARNING: no lineage found for tax_id $input_tax_id. Possible obsolete tax_id number.\n"; 
		next;
	}
	my $out_line = "$query{tax_id}\t";
	my $lineage = "";
	foreach my $key (sort {$b <=> $a} (keys(%Tax)))
	{
			$sql = "select name_txt from names where tax_id=$Tax{$key} and name_class like 'scientific name'";
			$sth = $dbh->prepare($sql)
			or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
			$sth->execute()
			or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");
			%query = %{$sth->fetchrow_hashref};
			
		# if anything was returned, check to be sure it's rank is approved
			my $current_name = $query{name_txt};				
			$sql = "select rank from nodes where tax_id=$Tax{$key}";
			$sth = $dbh->prepare($sql)
			or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
			$sth->execute()
			or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");
			%query = %{$sth->fetchrow_hashref};
			
			my $current_rank = $query{rank};
			if (exists $output_order{$current_rank} || $current_name =~/Viruses/)
			{
				#get rid of non-informative names	
					next if  $current_name =~ /unclassified .+/;
					$current_name =~ s/unclassified//;
					$current_name =~ s/ group//;
					$current_name =~ s/\(miscellaneous\)//;
					next if  $current_name =~ /candidate\s/;
					next if  $current_name =~ /candidate\s/;
					next if  $current_name =~ /\(class\)/;
					next if  $current_name =~ /subdivision/;
					next if  $current_name =~ /subgen/;
					next if  $current_name =~ /subg\./;
					next if  $current_name =~ /Family/;
					next if  $current_name =~ /cellular\sorganisms/;
					
					
				$lineage .= "$current_name;";
			}												
		}

	chop $lineage; # get rid of extra semicolon at end
	$|=1;
	if (defined $lineage && length $lineage >0)
	{
		print "$out_line"."$lineage\n";
	}
}

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

sub nodes(){
	my($tax_id,%Tax)=@_;
	$recursion_count++;
		
	$sql2 = "select parent_tax_id, rank from nodes where tax_id=$tax_id";
    $sth2 = $dbh->prepare($sql2)
            or dieStackTrace("Cannot prepare SELECT '$sql':\n<b>$DBI::errstr</b>");
    $sth2->execute()
        or dieStackTrace("Cannot execute SELECT '$sql':\n<b>$DBI::errstr</b>");
    %query2 = %{$sth2->fetchrow_hashref};
	$parental_tax_id=$query2{parent_tax_id};
	$Tax{$recursion_count}=$parental_tax_id;	
	
# end recursion when get to root level, or anything that has root as parent 
# ($recursion_max prevents runaway)

	if ($parental_tax_id == 1
		|| $parental_tax_id == 2   # Bacteria
		|| $parental_tax_id == 10239  #Viruses
		|| $parental_tax_id == 12884 # Viroids
		|| $parental_tax_id == 12908 # unclassified sequences
		|| $parental_tax_id == 28384  # other sequences
		|| $parental_tax_id == 131567  # cellular organisms
		|| $recursion_count > $recursion_max)  
	{
		$recursion_count = 0;
		$sth2->finish();
		return %Tax;
	}
	&nodes($parental_tax_id,%Tax);	
}

__END__

mysql> select * from names where name_class = "scientific name" and tax_id in (1,10239,12884,12908,28384,131567);
+--------+--------+------------------------+-------------+-----------------+
| id     | tax_id | name_txt               | unique_name | name_class      |
+--------+--------+------------------------+-------------+-----------------+
|      2 |      1 | root                   |             | scientific name |
|  36482 |  10239 | Viruses                |             | scientific name |
|  40841 |  12884 | Viroids                |             | scientific name |
|  40887 |  12908 | unclassified sequences |             | scientific name |
|  47466 |  28384 | other sequences        |             | scientific name |
| 241065 | 131567 | cellular organisms     |             | scientific name |
+--------+--------+------------------------+-------------+-----------------+
6 rows in set (0.00 sec)


select tax_id, parent_tax_id, rank from nodes where parent_tax_id = 1;
+--------+---------------+--------------+
| tax_id | parent_tax_id | rank         |
+--------+---------------+--------------+
|      1 |             1 | no rank      |
|  10239 |             1 | superkingdom |
|  12884 |             1 | superkingdom |
|  12908 |             1 | no rank      |
|  28384 |             1 | no rank      |
| 131567 |             1 | no rank      |
+--------+---------------+--------------+
6 rows in set (0.89 sec)

