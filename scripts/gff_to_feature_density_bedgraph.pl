#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use POSIX;

main();

sub main
{
	my $options = {};
	my $command_line = join(" ", ('gff_to_feature_density_bedgraph.pl', @ARGV));

	GetOptions($options, 'gff_file|g=s', 'tag|t=s', 'window|w=i', 'help|h');
	print STDERR keys(%$options) if( keys(%$options) != 3 );
	if( $options->{'help'} || keys(%$options) != 3 || !-f $options->{'gff_file'} || $options->{'tag'} =~ /\s/ || $options->{'window'} !~ /\d+/ ){
		print STDERR usage($command_line);
		exit;
	}

	my $windowSize = $options->{'window'};
	my %windows = (); #hash of windows
	my $fh = $options->{'gff_file'};
	open( GFF, "< $fh" ) or die "Cannot open GFF3 file: $fh\n";

	while( <GFF> ){
	
		chomp($_);
		next if( $_ =~ /^#.*$/ );
	
		#extract fields from each line of GFF3
		my ($chr, $source, $tag, $start, $end, $tr, $strand, $score, $notes) = split("\t", $_);
	

		if( $tag !~ $options->{'tag'} ){
			next;	
		}

		# Determine if feature rests on a window boundary:
		if( bin($start, $windowSize) != bin($start, $windowSize) ){
			#Feature is over a window boundary
			#Calculate number of bases in each window and increment totals

			#window boundary = bin# * $windowSize
			my $boundary_position = bin($start, $windowSize) * $windowSize;

			if( !defined($windows{$chr}{bin($start, $windowSize)}) ){ 
				$windows{$chr}{bin($start, $windowSize)} = 0;
			} elsif( !defined($windows{$chr}{bin($end, $windowSize)}) ){
				$windows{$chr}{bin($end, $windowSize)} = 0;
			}
			$windows{$chr}{bin($start, $windowSize)} += $boundary_position - $start + 1;
			$windows{$chr}{bin($end, $windowSize)} += $end - $boundary_position; #don't add 1 here (boundary_position is off by 1)
		} else {
			#Feature is within a single bin
			#Increment total bases covered for each window

			if( !defined($windows{$chr}{bin($start, $windowSize)}) ){ 
				$windows{$chr}{bin($start, $windowSize)} = 0;
			}
			$windows{$chr}{bin($start, $windowSize)} += $end - $start + 1;
		}

	}

	#Print the contents of $windows:
	foreach my $chr ( sort keys %windows ){
		foreach my $bin ( sort { $a <=> $b } keys %{$windows{$chr}} ){
			print STDOUT join("\t", $chr, ($bin * $windowSize - $windowSize + 1), $bin * $windowSize, $windows{$chr}{$bin}/$windowSize) . "\n"; 
		}
	}
}

#given an integer and windowSize, calculate the binNumber
sub bin
{

	my ($pos, $windowSize) = @_;

	return floor($pos/$windowSize) + 1;

}

sub usage
{
	my ($command) = @_;
  
  my $usage .= <<__HERE__;

  Required Parameters:

    -g, --gff-file : filename for the gff3 containing features
    -t, --tag : the tag for the feature to use when calculating density
    -w, --window: the size of windows to calculate feature density

  Other Parameters:

    -h, --help : display this usage information

  Your command:

    $command

__HERE__

return $usage;
}
