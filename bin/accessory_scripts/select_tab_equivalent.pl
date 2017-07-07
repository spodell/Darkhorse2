#!/usr/bin/env perl 
# select_tab_equivalent.pl 
# Sheila Podell 
# October 21, 2005

#	Selects all lines from an input tabfile having same primary key value and
#	"equivalent" values (within tolerance) in a second specified field.

#	"Equivalent" is defined as being within a tolerance range where difference
#	is considered insignificant

# 	input parameters:
#		infile name
#			string (valid path, readable, non-zero size)
#		primary key field number (default = 0, first field)
#			positive integer < total number of fields
#		selection field number (default = 1, second field)
#			positive integer < total umber of fields AND != primary key
#		minimum/maximum (default maximum)
#			string "min" or "max"
#		percentage/absolute = tolerance as pct or abs value (default = pct)
#			string "pct" or "abs"
#		tolerance = difference considered equivalent
#			positive float between 0 and 1

#	output - new file containing a subset of lines from input file
# MODIFIED 4/12/06 to send output to STDOUT
# MODIFIED 5/28/06 to use <= instead of < as filtration criterion
# this fixes bug where zero tolerance didn't work (now returns tophit)

use warnings; 
use strict;
use Getopt::Long;

# Default values
	my $in_filename; # = "test.tab";
	my $key_field = 0;
	my $select_field;  # = 1;
	my $min_max = "max";
	my $pct_abs = "pct";
	my $tolerance = 0.1;
	my $limit = 100;
	my $USAGE;
	my $message = qq(   Usage: $0 -i infile_name  -s select_field <additional options>;
	
  Required parameters
  	-i  <input file name>
  	-s	<selection field number>  [default = 1] 
  
  Optional parameters:
	-k	<primary key field number> [default = 0]
	-m	<min> or <max> [default max] - select minimum or maximum values
	-p	<pct> or <abs> [default percent] - select percentage or absolute tolerance
	-t	<tolerance> [default 0.1] - floating point number
	-l  <limit number of selections per key> [default = 100]
	
	Description:
		Selects all lines from input tabfile having same primary key value with
		"equivalent" values (within tolerance) in a second specified field.
	)."\n\n";

# get command line arguments, check validity
	GetOptions(     "i=s"	=> \$in_filename,		# infile name
					"k=i"	=> \$key_field,		# primary key field number
					"s=i"	=> \$select_field,	# selection field number
					"m=s"	=> \$min_max,		# use minimum or maximum for selection
					"p=s"	=> \$pct_abs,		# use percent or absolute values for tolerance
					"t=f"	=> \$tolerance,		# tolerance (difference considered insignificant)
					"l=i"	=> \$limit,			# limit number of selections per key
        );        
    
	if($USAGE || !$in_filename || !$select_field) 
    {
		print $message;
		exit(0);
	}
	unless (-s $in_filename)
	{
		print "cannot open infile $in_filename\n";		
		exit (0);
	}	
	unless ($key_field > -1 && $key_field < 257)
	{
		print "selection field has invalid column range (must be 0-256)\n";		
		exit (0);
	}
	unless ($select_field > -1 && $select_field < 257)
	{
		print "selection field has invalid column range (must be 0-256)\n";
		exit (0);
	}
	unless ($min_max eq "min" || $min_max eq "max")
	{
		print " \n  field -m must be specified as \"min\" or \"max\"\n\n";
		print "$message";
		exit (0);
	}
	unless ($pct_abs eq "pct" || $pct_abs eq "abs")
	{
		print "field -p must be specified as \"pct\" or \"abs\"\n";
		print "for tolerance by percentage or absolute value\n";
		exit (0);
	}
	unless ($limit > -1)
	{
		print "limit must be a positive integer.\n";		
		exit (0);
	}

# open filehandles
	my $outfile = substr ($in_filename, 0, 10);
	$outfile .= ".select";
	open (INPUT,"$in_filename") or die "can't open input file $in_filename\n$!\n";	
	#open (OUTPUT, ">$outfile") or die "can't open outfile $outfile\n $!\n";
	
# group the lines into arrays, one for each different primary key
	my %keylist = ();	# list of unique keys, each associated with an array of lines
	my %keylist2 = ();	# list of unique keys, each associated with a data structure
	my @ordered_keylist = ();
	my %headers = ();
	open (INPUT,"$in_filename") or die "can't open input file $in_filename \n$!\n";	
	while (<INPUT>) 
	{
		next if ($_ =~ /^\s/);
		chomp;
		my @params = split "\t", $_;		
		my $key = $params[$key_field];	
	
		unless (exists $keylist{$key}) 
		{
		# keep track of order in original file
			push @ordered_keylist, $key;
		# create a new array to hold all lines containing same key
			my @array = (); 
  			$keylist{$key} = \@array;	
		} 	
	
	# add current line to the list of lines associated with that key 
		push @{$keylist{$key}}, $_;
		
	}
	close INPUT;
	
#  sort keylist, subsort lines by selection field, select those that meet criteria
	foreach (@ordered_keylist)
	{
		next if ($_ =~ /query/); #skip header line
		my @sorted_lines = &sort_by_select_field(\@{$keylist{$_}}); #puts max or min first
		my @initial_params = split "\t", @{$keylist{$_}}[0];
	
	# deal with corrupted data
		unless (defined $initial_params[$select_field] && $initial_params[$select_field] =~ /[\de\.-]+$/)
		{			
			next if ($. == 0); #header
			print STDERR "no quantitative data in select field for the following line:\n";
			print STDERR "  @{$keylist{$_}}[0]\n";
			next;	
		}		
		
		my $initial_select = 0 + $initial_params[$select_field]; # use 0 + to ensure scalar		
		
		my $num_selected = 0;
		foreach my $line (@sorted_lines)
		{
			my @current_params = split "\t", $line;
			my $current_select = 0 + $current_params[$select_field];			
			my $current_diff = abs($current_select - $initial_select);
		
		#	hack to deal with zero values for selection column (e.g. blast e-values)
			my $pct_diff = 0;
			if ($initial_select == 0)
			{
				$initial_select = 0 + 1e-250;
			}			
			
			$pct_diff = $current_diff/$initial_select;
			
		# keep only lines that fall within tolerance range			
			if ($pct_abs eq "pct" && $pct_diff <= $tolerance)
			{
				print  "$line\n";
				$num_selected++;
			}
			if ($pct_abs eq "abs" && $current_diff <= $tolerance)
			{
				print  "$line\n";
				$num_selected++;
			}
			
			last if ($num_selected > $limit -1);
		}
	}
	
# user output
	close INPUT;	


#####################
# SUBROUTINES
#####################

sub new {
  my ($className, $param) = @_;
  my $self = {};
  bless $self, $className;
  my @properties = @$param;
  $self->{name} = $properties[0]; 
  $self->{lines} = $properties[1]; 
  $self->{select_value} = ""; 
  $self->{orig_order} = "";

  return($self)
}


sub sort_by_select_field
{
	my ($array_ref) = @_;
	my @lines = @$array_ref;
	my %line_sorter = ();
	foreach my $line (@lines)
	{
		my @tmp =split "\t", $line;
		
	# deal with corrupted data
		unless (defined $tmp[$select_field] && $tmp[$select_field] =~ /[\de\.-]+$/)
		{
			print STDERR "no quantitative data in select field for the following line:\n";
			print STDERR "  @{$keylist{$_}}[0]\n";
			next;	
		}
		$line_sorter{$line} = $tmp[$select_field];			
	}
		
	my @sorted_array = ();
	if ($min_max eq "min")
	{
		@sorted_array = sort {$line_sorter{$a} <=> $line_sorter{$b}} keys %line_sorter;
	}
	else
	{
		@sorted_array = sort {$line_sorter{$b} <=> $line_sorter{$a}} keys %line_sorter;
	}
	
	
	return @sorted_array;
}

		
__END__
