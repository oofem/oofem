#!/usr/bin/perl
#
# set RESTRICTED on itself (this file will not be included in public release)
# <restricted_file> 
#
# utility to exclude from oofem source the
# sections and files, which are not supposed
# to be in (public) release. 
#
# WARNING:
#
# The original file is overridden !!!!
#
# each file is scanned for the following keywords:
#
# 1) <restricted_section> and </restricted_section>
#    The section contained within these keywords will be removed
#
# 2) <restricted_file> 
#    The whole processed file will be excluded (removed)
#
#
# written by Borek Patzak (c) 2001
#
$proceed = 0;
$exclude = 1;
$skip_file = 0;
$parse_mode = $proceed;

if (@ARGV != 1)
{
		printf ("Usage: release_filter.pl filename\n");
} 



$fileName = $ARGV[0];
open (INPUT, $fileName)
		or die "Can't open $fileName for reading: $!";

$tmpfilename = "tmpfile";
open (TMPFILE, ">$tmpfilename")
		or die "Can't open $tmpfilename for writing: $!";

#parse the file and serach for keywords
while (($line = <INPUT>) && ($skip_file == 0))
{
#		chop $line;
		&parseLine ($line);
}

close (TMPFILE);

if ($skip_file == 0) 
{
		if ($parse_mode == $exclude) {
				die ("error: exclude_section missing end\n");
		}
		rename ($tmpfilename, $fileName)
				or die "Can't rename $tmpfilename to $fileName: $!";
} else {
		print "exclude_file detected: removing $fileName\n";
		unlink $fileName;
}


sub parseLine {
		local ($line) = @_;
		local ($strBefore, $strAfter);


#		print "GET: $line";

		if ($skip_file) 
		{
				return;
		}

		if ($line =~ /(.*)<(restricted_section|RESTRICTED_SECTION)>(.*)/) {
				$strBefore = $1;
				$strAfter  = $3;

				&parseLine ($strBefore);
				if ($parse_mode == $exclude) {
						die "error: nested exclude_sections detected: $fileName\n";
				}
				$parse_mode = $exclude;
				&parseLine ($strAfter);
		} elsif ($line =~ /(.*)<\/(restricted_section|RESTRICTED_SECTION)>(.*)/) {
				if ($parse_mode == $proceed) {
						die "error: end of non-opened exclude_sections detected\n";
				}
				$strBefore = $1;
				$strAfter  = $3;

				$parse_mode = $proceed;
				&parseLine ($strAfter);
				print TMPFILE "\n"; 
		} elsif ($line =~ /(.*)<(restricted_file|RESTRICTED_FILE)>(.*)/) {
#				print "exclude_file detected\n";
				$skip_file = 1;
				return;
		} else {

				if ($parse_mode == $proceed)
				{
						print TMPFILE $line;
				}
		}
}
		

