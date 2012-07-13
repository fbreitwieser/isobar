#!/usr/bin/perl
# Creation date : 2011-03-28
# Last modified : Wed 27 Jun 2012 02:14:32 PM CEST

# Module        : psx2tab2.pl
# Purpose       : 
# Usage         : 
# Licence       : Copyright (c) 2010 Florian Breitwieser, Ce-M-M-
# Contact       : Florian Breitwieser <fbreitwieser@cemm.oeaw.ac.at>

use strict;
use warnings;
use DBI;

my ($protID, $peptide, $modif, $charge, $theo_mass, $exp_mass, 
    $parent_intens, $start_pos, $rtime, @search_engine, @score,$spectrum);

my $numb = "[0-9]+\.?[0-9]*(?:e[+-][0-9]*)?";

print "accession\tpeptide\tmodif\tcharge\ttheo.mass\texp.mass\tparent.intens".
      "\tstart.pos\tretention.time\tsearch.engine\tscore\tspectrum\n";

foreach my $psxFile (@ARGV) {
eval{
  my $XML;
  my $is_stdin = 0;
  if (defined $psxFile){
    open $XML, "<", $psxFile or die $!;
  } else {
    $XML = *STDIN;
    $is_stdin++;
  }
  while (<$XML>){
    if (/<idi:proteinId>(.*)<\/idi:proteinId>/)    { $protID = $1; 
    } elsif (/<idi:sequence>(.*)<\/idi:sequence>/) { $peptide = $1;
    } elsif (/<idi:theoMass>(.*)<\/idi:theoMass>/) { $theo_mass = $1; 
    } elsif (/<idi:charge>(.*)<\/idi:charge>/)     { $charge = $1; 
    } elsif (/<idi:startPos>(.*)<\/idi:startPos>/) { $start_pos = $1; 
    } elsif (/<idi:modif>(.*)<\/idi:modif>/)       { $modif = $1; 
    } elsif (/<idi:retentionTime>(.*)<\/idi:retentionTime>/) { $rtime = $1; 
    } elsif (/<ple:PeptideDescr>.!.CDATA.([^"]+)..><\/ple:PeptideDescr>/) {
        $spectrum = $1;
    } elsif (/<ple:ParentMass><!.CDATA.($numb) ($numb) [0-9\.,]+..><\/ple:ParentMass>/) {
        $exp_mass = $1;
        $parent_intens = $2;
    } elsif (/<idi:peptScore engine="([^"]+)".*>(.*)<\/idi:peptScore>/) {
        push @search_engine, $1;
        push @score,$2;
    }
    elsif (/<\/idi:OneIdentification>/){
        die "not defined" if (!defined $spectrum);
        print join ("\t",$protID, $peptide, $modif, $charge, 
                    $theo_mass, $exp_mass, $parent_intens, 
                    $start_pos, $rtime, join("|",@search_engine), join("|",@score), $spectrum)."\n";
    
        undef $peptide;
        undef $modif;
        undef $charge;
        undef $theo_mass;
        undef $exp_mass;
        undef $parent_intens;
        undef $start_pos;
        undef $rtime;
        undef @search_engine;
        undef @score;
        undef $spectrum;

    }
  }
  close($XML) unless $is_stdin;
};

if ($@){
  die("Encountered a problem: $@");
}
}

