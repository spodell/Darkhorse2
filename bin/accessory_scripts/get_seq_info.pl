#!/usr/bin/perl
# get_seq_info.pl
# Sheila Podell 
# May 21, 2006

#	takes in a sequence file in fasta format 
#	outputs to STDOUT a tab-delimited table:
#	sequence_id	  seq_length   annotation

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
	my $flag = 0;
	my @seq_list = ();  
	my $current = "";
	
# get sequence, id, annotation
	while (my $line = <INFILE>)
	{
		next if $line =~ /^\s+$/;
		chomp $line;		
		if ($line =~ />(.+)/) 
		{								
			my @tmp = split " ", $line;
			$current = new_seq("record",\@tmp);
			push (@seq_list, $current);
			$flag = 1;
		}													
		else
		{			
			$current->{sequence} .= uc"$line";
		}	
	}
	close INFILE;	
	foreach my $next (@seq_list)
	{
		my $seq_listtring = $next->{sequence};
		$next->{seq_length} = length ($seq_listtring);
		print  "$next->{id}\t$next->{seq_length}\t$next->{annotation}\n";
	}
	
###################################################
# SUBROUTINES
###################################################
sub new_seq {
  my ($className, $param) = @_;
  my $self = {};
  bless $self, $className;
  my @properties = @$param;
  my $header = join " ", @properties;
  $self->{header} = $header;  
  $properties[0] =~ s/\>//;
  $self->{id} = shift @properties;   
  $self->{annotation} = join " ", @properties;
  $self->{sequence} ="";
  $self->{seq_length} ="";  
  return($self)
}