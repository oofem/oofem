#!/usr/bin/perl

$npar=$#ARGV;
@formats=@ARGV;

if ($npar == -1) {
		&printHelp();
}

while ($line = <STDIN>)
{
		@c = split(' ',$line);
		foreach $expr (@formats) {
#				print $expr,"=",eval($expr)," ";
				print eval($expr)," ";
		}
		print "\n";

}

sub printHelp {
		print "\nusage: summer 'expr1' 'expr2' ... 'exprN'\n\n";
		print "where 'expr' is any valid math perl expression (use the ''\n";
		print "to prevent from shell substitution).\n";
		print "To access particular columns of input file (stdin) use c array\n";
		print '($c[0] for first column value, $c[1] for second, etc.)',"\n";
		print "The output will contain evaluated expressions separated by space.\n";
		print "\n\nAuthor: Borek Patzak\n";
		exit;
}
