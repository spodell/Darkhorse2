#!/usr/bin/env perl
# getseq_multiple.pl
# Sheila Podell 
# December 8, 2004

# inputs:
#	a file containing multiple fasta format sequences
# 	a list of accession numbers to retrieve (newline separated). 
# outputs selected records to STDOUT
#
# suitable for really large (gigabyte) fasta source files
# and large numbers of retrieval sequences

use warnings; 
use strict;

# get filenames
	unless (scalar @ARGV ==2)
	{
		print STDERR "Usage: $0 id_list fasta_filename\n";
		exit(0);
	}

# get selection keys
	my %selections = ();
	my $selectname = $ARGV[0];
	open (SELECT, "<$selectname") or die "couldn't open $selectname, $!\n";
	while (<SELECT>)
	{
		chomp;
		$selections{$_} = 0;
	}	
	close SELECT;
	my $num_requests = scalar (keys %selections);

# go thru infile (once only), select all records that match list 
	my $infilename = $ARGV[1];
	open (INFILE, "<$infilename") or die "couldn't open $infilename, $!\n";
	my $select = 0;		# 0=false, 1= true
	my $count_hits = 0;
	
	while (my $line = <INFILE>)
	{
		if ($line =~ /^\>(.+)/)
		{
			last if ($count_hits == $num_requests);
			$select = &check_name($1);
			if ($select > 0)
			{
				print $line;
				$count_hits++;
			}
			next;
		}				
		if ($select > 0)
		{
			print $line;	
		}		
	}

	
	close INFILE;
		
#################################
# SUBROUTINES
#################################

sub check_name
{
	my ($header) = @_;
	my $result = 0;
	chomp $header;
	my @seq_words = split (" ", $header);
	my $seq_id = $seq_words[0];
		
	if (exists $selections{$seq_id})
	{
		$result = 1;
	}
	return $result;		
}

__END__


