#!/usr/bin/env perl
# get_nr_seq_info_v2.pl
# Sheila Podell 
# October 2, 2016

#	takes in a sequence file in fasta format 
#	outputs to STDOUT a tab-delimited table:
#	sequence_id	  seq_length   annotation

#  minimizes memory usage for really huge files
# limits annotation to 80 characters, and removes weird punctuation

# February 6, 2012 - fixed bug that failed to print last fasta sequence 
# if end of file was a blank line.

# October 2, 2016 - removed dependence on gi_num

use warnings; 
use strict;

# get in & out file handles
	unless (@ARGV == 1)
	{
		print "usage: $0 fasta_filename\n";
		exit (0);
	}	
	my $infilename = $ARGV[0];
	open (INFILE, "<$infilename") or die "couldn't open $infilename, $!\n";	 
	my $current = "";
	
# get sequence, id, annotation
	while (my $line = <INFILE>)
	{
		chomp $line;
		if (eof) 
		{					
			unless ($line =~ /^\s+$/)
			{
				$current->{sequence} .= uc"$line";
			}
			
			my $id = $current->{id};
				my $seq_length = length $current->{sequence};
				my $annotation = $current->{annotation};
				my @printarray = (
					$id,
					$annotation,
					$seq_length
					);
				my $printstring = join "\t", @printarray;
				print "$printstring\n";
				last;
		}
		
		next if $line =~ /^\s+$/;		
		if ($line =~ />(.+)/ || eof) # this prints last fasta at end of file without blank line
		{								
		# print previous object unless it hasn't been created yet
			unless ($. < 2)
			{				
				my $id = $current->{id};
				my $seq_length = length $current->{sequence};
				my $annotation = $current->{annotation};
				my @printarray = (
					$id,
					$annotation,
					$seq_length
					);
				my $printstring = join "\t", @printarray;
				print "$printstring\n";				
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
	
	
###################################################
# SUBROUTINES
###################################################
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
 
 # limit annotation to a reasonable size	
	$self->{annotation} = substr($self->{annotation}, 0, 80); 
  
  $self->{sequence} ="";
  $self->{seq_length} =""; 
  
  return($self)
}