#!/usr/bin/perl

#
# utility to find last solution step time from oofem output file
# input : oofem output file name
# output: last solution step time written to stdout
# written by Borek Patzak (c) 2000
#

if (@ARGV != 1)
{
		printf ("Usage: findlaststeptime oofemOutFile \n");
} 

$oofemOutFile = $ARGV[0];

open (OOFEMOUT, $oofemOutFile);

#parse output file and serch solution step header
$lastStepTime = 0;

while ($line = <OOFEMOUT>) 
{
		chop $line;
		if ($line =~ /^O/) #for fast record elimination
		{
				if ($line =~ /^Output for time\s+(\d+.\d+e[+-]\d+).*solution step number\s+(\d+)/) 
				{
						$lastStepTime = $1;
						$lastStepNumber = $2;
				}
		}
}

#print  int $lastStepTime;
print  int $lastStepNumber;
