#!/usr/bin/perl
# Creation date : 2012-03-18
# Last modified : Sun 18 Mar 2012 11:53:47 AM CET

# Module        : tab2psx.pl
# Purpose       : 
# Usage         : 
# Licence       : Copyright (c) 2010 Florian Breitwieser, Ce-M-M-
# Contact       : Florian Breitwieser <fbreitwieser@cemm.oeaw.ac.at>

use strict;
use warnings;
use DBI;
use Data::Dumper;

my $idfile = shift @ARGV;

my $numb = "[0-9]+\.?[0-9]*(?:e[+-][0-9]*)?";

my $ID;
open $ID, "<", $idfile or die $!;
 
my $header = <$ID>; chomp $header;
my $i=0;
my %header = map { $_ => $i++ } (split(/\t/,$header));
print STDERR Dumper \%header;
#die;

my $date = `date +%Y-%m-%d`; chomp $date;
my $time = `date +%H:%M:%S`; chomp $time;

print <<eoi;
<?xml version="1.0" encoding="ISO-8859-1"?>
  <idi:ProtSpectraIdentifications  xmlns:idi="namespace/ProtSpectra.html">
    <idi:OneSample>
      <idi:header>
        <idi:instrument>n/a</idi:instrument>
        <idi:spectrumType>msms</idi:spectrumType>
        <idi:date>$date</idi:date>
        <idi:time>$time</idi:time>
        <idi:searchEngines>
          <idi:oneEngine />
          </idi:oneEngine>
        </idi:searchEngines>
        <ple:ItemOrder xmlns:ple="namespace/PeakListExport.html">
          <ple:item type="mass"/>
          <ple:item type="intensity"/>
          <ple:item type="charge"/>
        </ple:ItemOrder>
      </idi:header>
    <idi:Identifications>
eoi

my $oldprotein;
while (<$ID>) {
  chomp;
  my @data = split(/\t/);
  my $protein = $data[$header{'accession'}];
  if (!defined $oldprotein || $protein ne $oldprotein) {
    if (defined $oldprotein) {
      print "      </idi:OneProtein>\n";
    }
    print "      <idi:OneProtein>\n";
    print "      <idi:proteinId>$protein</idi:proteinId>\n";
    $oldprotein = $protein;
  }
  
  print
"      <idi:OneIdentification>
        <idi:answer>
          <idi:sequence>",$data[$header{'peptide'}],"</idi:sequence>
          <idi:modif>",$data[$header{'modif'}],"</idi:modif>
          <idi:theoMass>",$data[$header{'theo.mass'}],"</idi:theoMass>
          <idi:charge>",$data[$header{'charge'}],"</idi:charge>
          <idi:startPos>",$data[$header{'start.pos'}],"</idi:startPos>
          <idi:retentionTime>",$data[$header{'retention.time'}],"</idi:retentionTime>
        </idi:answer>
        <idi:source>
          <idi:file>$idfile</idi:file>\n";
  if (defined $header{'score'} && defined $header{'search.engine'}) {
    my @se = split(/\|/,$data[$header{'search.engine'}]);
    my @score = split(/\|/,$data[$header{'score'}]);
    for (my $j=0;$j<=$#se;++$j) {
      print "          <idi:peptScore engine=\"",$se[$j],"\">",$score[$j],"</idi:peptScore>\n";
    }
  }
  print 
"        </idi:source>
        <ple:peptide spectrumKey=\"",$data[$header{'spectrum'}],"\" xmlns:ple=\"namespace/PeakListExport.html\">
        <ple:PeptideDescr><![CDATA[",$data[$header{'spectrum'}],"]]></ple:PeptideDescr>
        <ple:ParentMass><![CDATA[",$data[$header{'exp.mass'}]," ",$data[$header{'parent.intens'}]," ",$data[$header{'charge'}],"]]></ple:ParentMass>\n";
  print       
"        <ple:peaks><![CDATA[
]]></ple:peaks>
        </ple:peptide>
      </idi:OneIdentification>\n";
}
print <<eoi;
    </idi:Identifications>
  </idi:OneSample>
</idi:ProtSpectraIdentifications>
eoi

close $ID;
