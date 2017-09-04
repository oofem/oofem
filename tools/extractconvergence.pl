#!/usr/bin/perl


#
# utility to extract from oofem output the convergence 
# data for each step.
# for each step the file c$i.dat is created, where $i is step number
# for each iteration, the follofing line is printedL
# #iter #loadLevel #forceError #dispError
#
# written by Borek Patzak (c) 2000
#

#read from stdin
open (INPUT, '-');
$closeflag = 0;

print ("found steps: ");

while ($line = <INPUT>)
{
		chop $line;
		if ($line =~ /Solving \[step number\s+(\d+)\]/)
		{
				if ($closeflag == 1) 
				{
						close (OUTPUT);
				}
				open (OUTPUT, ">c$1.dat");
				print "$1 ";
				$closeflag = 1;
		}
		if ($line =~ /(\d+)\s+(\d+.\d+e[+-]\d+)\s+(\d+.\d+e[+-]\d+)\s+(\d+.\d+e[+-]\d+)/)
		{
				print OUTPUT "$line\n";
		}
}

print ("\nThanks for using extractconvergence.pl\n");
