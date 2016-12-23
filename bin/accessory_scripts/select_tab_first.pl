#!/usr/bin/perl 
# select_tab_first.pl 
# Sheila Podell 
# August 2, 2005

#	Takes in a file in tab format 
#	selects first instance of each category in column 0
# (assumes file has been pre-sorted by some other criteria
# in cases of multiple occurrences for col 0)

# results to STDOUT

use warnings; 
use strict;

my %keywords = ();
my $num_tablines = 0;
my $num_hits = 0;

# Get filenames & open filehandles 
	unless (scalar @ARGV == 1)
	{
		print "Usage: $0 tab_filename\n";
		exit (0);
	}
	
	my $in_filename = $ARGV[0];
	my $out_filename = substr ($in_filename, 0, 10);
	$out_filename .= ".firsts";
	open (INPUT,"$in_filename") or die "can't open input file\n$!\n";			
	
# Check for keys	
	while (<INPUT>) 
	{
		next if ($_ =~ /^\s/);
		chomp;
		$num_tablines++;
		my @tmp = split "\t", $_;
		unless (exists $keywords{$tmp[0]})
		{
			$keywords{$tmp[0]}++;
			$num_hits++;
			print "$_\n";
		}				
	}	
	
#user output
	my $num_kw = scalar (keys %keywords);
	
	close INPUT;	
	print "\n";
		
	
__END__
