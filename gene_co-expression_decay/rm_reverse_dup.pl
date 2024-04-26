#!/usr/bin/bash perl
use strict;
use warnings;

my %store;
my $header = 'col	row	value.x	value.y';
while (<>){
	chomp;
	s/^[0-9]+\s+//;
	my ($A, $B) = (split)[0,1];
	my $AB = "${A}_$B";
	my $BA = "${B}_$A";
	next if /^col\s+row\s+/;
	$store{$AB} = $_ unless exists $store{$AB} or exists $store{$BA};
}

print "$header\n";
foreach my $key (keys %store){
	print "$store{$key}\n";
	}

