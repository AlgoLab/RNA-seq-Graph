#!/usr/bin/perl
use bigint;

my $r= "";
while (<>) {
	chomp;
	my $s= $_;
	foreach my $c (split //, $s) {
		if ($c eq "A") {
			$r .= "00";
		} elsif ($c eq "C") {
			$r .= "01";
		} elsif ($c eq "G") {
			$r .= "10";
		} else {
			$r .= "11";
		}
	}
	my $i= 0;
	foreach my $c (split //, $r) {
		$i *= 2;
		if ($c eq "1") {
			$i ++;
		}
	}
	print $i,"\n";
}
