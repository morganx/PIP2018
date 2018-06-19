#! /usr/bin/perl-w
use strict;
use Data::Dumper;

my %kegg;

get_kegg();
#print Dumper %kegg;
while(<>)
{
chomp;
my @l = split/\t/;
my $a = $l[2];
$a = $a." $kegg{$a} ";
$_=~s/$l[2]/$a/;
print $_, "\n";
}

sub get_kegg
{
open(FILE, "~/PIP2018/src/kegg.txt") or die "could not get kegg annots";
while(<FILE>)
	{
	chomp;
	my($id, $desc) = split(/\t/);
	$kegg{$id} = $desc;
	}
close(FILE);	
}