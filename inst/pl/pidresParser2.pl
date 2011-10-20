#!/usr/bin/perl

# Mass spectrometry Perl program for extracting correct peptide matches from Phenyx pidres.xml files

# Copyright (C) 2008, 2009, 2010 Jacques Colinge

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.

# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# Contact:
#  Jacques Colinge
#  CeMM
#  Lazarettgasse 19/3
#  A-1090 Vienna, Austria
#  www.cemm.at

=head1 NAME

pidresParser2.pl - Extraction of reliable protein/peptide/spectrum matches from Phenyx pidres.xml files

=head1 SYNOPSIS

pidresParser2.pl --defaultset=SETNAME [options] pidres.xml files

=head1 OPTIONS

Use pidresParser2.pl -h

=head1 DESCRIPTION

The script parses one or several Phenyx pidres.xml files to extract reliable peptide/spectrum matches and
outputs them in the .protSpectra.xml format. The pidres.xml file(s) can be compressed (gzipped) files.

The selection of the peptide assignments is performed based on several thresholds applied to identifications
found in the pidres.xml file(s):

=over 4

=item -help -h

=item -verbose

=item -parseall

Parses all the files, otherwise BSA, BLK, and ACN runs are skipped.

=item -notonlybestpept

Not only keeps only the best peptide for each spectrum but all.

=item -distinctprot

Only outputs proteins with at least one specific peptide.

=item -nofilelevel

Global instead of file-based check of the selection criteria. Useful for considering several chromatographic/gel fractions as a single sample.

=item --definitionfile=fname

The name of an XML file containing amino acid and modification mass definitions. By default, the parser reads a file names parsersConfig.XML in the directory where it is itself located.

=item --proteinforce=AC-list

Let specific proteins pass through any threshold (but the basicscore) for QC or increased sensitivity, AC-list must be comma-separated. One use is to follow specific peptides in AP tags.

=item --minbigred=int

Number of specific peptides to be a big red (see Mascot documentation); only proteins having a big red peptide are selected.

=item --minscore=float

Minimum peptide z-score.

=item --maxpvalue=float

Maximum peptide P-value.

=item --minprotscore=float

Minimum protein score.

=item --minnumpept=int

Minimum number of distinct peptides per protein.

=item --savescore=float (single peptide hits)

Minimum peptide save z-score to bypass the minimum number of peptides required to be above minscore (allows single peptide hits with a higher threshold). Default is +infinity.

=item --basicmodiflist=coma separated list

List of basic modifications that should just be processed with the normal thresholds. Default is 'Cys_CAM,Oxidation_M'.

=item --specialmodifscore=float

Minimum z-score for the modified peptides whose modifications do not all belong to the basicmodiflist (applies to both validation and output). Default is +infinity, i.e. such peptides are filtered out.

=item --minseqcovsph=float

Minimum sequence coverage for single peptide hits

=item --minlen=int

Minimum peptide sequence length.

=item --defaultset=['cid','cid_sensitive','hcd','qtof','qtof_sensitive']

Set the parameters to defaults according to a set of predifined values. The default values set not using this parameter are more conservative. The command line changes overwrite the defaults. Hence it is possible to take a default set and to change a specific parameter through the command line.

=item -lightxml

Flag to avoid output mass lists (much small XML)

=back

To be selected a peptide must have a z-score larger than the minimum peptide score, a p-value smaller than the maximum p-value. The protein is finally selected if it has minnumpept peptides above minscore at least and a protein score better than minprotscore.

Alternatively, it is possible to specify savescore and then the condition on the number of peptides with z-score better than minscore is substitued by the following one: sum of the peptide z-scores better than minscore must be above savescore and a minimum sequence coverage of minseqcovsph must be reached. Sequence coverage is determined using all the peptides > outputscore. This alternative validation is mostly used for single peptide hits (SPHs).

In every case, as soon as a protein is selected, all the peptides above outputscore (usually < minscore) are reported.

During the parsing of the file, each spectrum is associated with the peptide that gives the best match, i.e. all multiple interpretations of a spectrum are lost in favor of the best one, unless -notonlybestpept flag is set. Moreover, all peptides with score less than the basic score (typically 5) are not read.

=head1 EXAMPLE

./pidresParser2.pl example.pidres.xml > test.protSpectra.xml

=head1 AUTHOR

Jacques Colinge

=cut

use strict;
use warnings;
no warnings 'uninitialized';
use Config::IniFiles;
use FindBin qw/$Bin/;
use File::Basename;
use Getopt::Long;
use Pod::Usage;

my $phenyxVersion = '2.5.14';
my $internSequence = 0;

my $massType = 'monoisotopic'; # Hardcoded at CeMM
my ($database, $dbRelease, $dbReleaseDate, $fixedModifs, $variableModifs, $nmc, $enzymeName, $cleaveMode);
my ($searchFragmentTol, $searchFragmentTolUnit, $searchParentTol, $searchParentTolUnit);

# To compliant with older cemmpidres where the scoring function is not repeated in the output. To be removed asap
$searchFragmentTol = '0.5';
$searchFragmentTolUnit = 'Da';

my ($help, $verbose);
my $instrument = 'n/a';
my $definitionFile = $Bin."/parsersConfig.xml";
my $configFile = $Bin."/pidresParser.ini";
my $qrValidAASeq = qr/^[ACDEFGHIJKLMNOPQRSTUVWY]+$/;
my ($notOnlyBestPept, $distinctProt, $noFileLevel, $parseAll, $defaultSet, $lightXML, $proteinForce);

my $config = Config::IniFiles->new(-file=>$configFile) 
    or die ("Could not intialize configuration file: $!");

# Checks whether a profile is given
for (my $i = 0; $i < @ARGV; $i++){
  if (index($ARGV[$i], '-defaultset=') != -1){
    $ARGV[$i] =~ /\-defaultset=(.+)/;
    $defaultSet = $1;
    splice(@ARGV, $i, 1);
    last;
  }
}

if (!defined($defaultSet) || !$config->SectionExists(uc($defaultSet))) {
  print STDERR "Use on of the available parameter sets using --defaultSet=SETNAME:\n".
      "\t",join("\n\t",$config->Sections())."\n";
  exit(1);
}

my $minScore          = get_val($defaultSet,'minScore');
my $outputScore       = get_val($defaultSet,'outputScore');
my $maxPvalue         = get_val($defaultSet,'maxPvalue');
my $defaultSaveScore  = get_val($defaultSet,'defaultSaveScore');
my $saveScore         = get_val($defaultSet,'saveScore');
my $specialModifScore = get_val($defaultSet,'specialModifScore');
my $basicModifList    = get_val($defaultSet,'basicModifList');
my $minSeqCovSPH      = get_val($defaultSet,'minSeqCovSPH');
my $minProtScore      = get_val($defaultSet,'minProtScore');
my $minNumPept        = get_val($defaultSet,'minNumPept');
my $minLen            = get_val($defaultSet,'minLen');
my $minBigRed         = get_val($defaultSet,'minBigRed');

my $result = 
  GetOptions('help' => \$help,
             'h' => \$help,
		'parseall' => \$parseAll,
		'notonlybestpept' => \$notOnlyBestPept,
		'distinctprot' => \$distinctProt,
		'lightxml' => \$lightXML,
		'definitionfile=s' => \$definitionFile,
		'nofilelevel' => \$noFileLevel,
		'proteinforce=s' => \$proteinForce,
		'savescore=f' => \$saveScore,
		'minseqcovsph=f' => \$minSeqCovSPH,
		'minbigred=i' => \$minBigRed,
		'outputscore=f' => \$outputScore,
                'specialmodifscore=f' => \$specialModifScore,
                'basicmodiflist=s' => \$basicModifList,
		'minscore=f' => \$minScore,
		'minnumpept=i' => \$minNumPept,
		'minprotscore=f' => \$minProtScore,
		'maxpvalue=f' => \$maxPvalue,
		'minlen=i' => \$minLen,
		'instrument=s' => \$instrument,
        'verbose' => \$verbose);

pod2usage(-verbose=>2, -exitval=>2) if (!$result);
pod2usage(get_params()) if (defined($help) && defined($defaultSet));
pod2usage(-verbose=>2, -exitval=>2) if (defined($help));
pod2usage(".pidres.xml files must be provided on command line!") if (!defined @ARGV || skalar(@ARGV) == 0);

if (defined($outputScore)) {
  pod2usage("outputScore (for hits of validated proteins) may not be higher than minScore (for validation)!") if ($outputScore > $minScore);
}

$outputScore = $minScore if (!defined($outputScore));

my %proteinForce;
if (defined($proteinForce)){
  foreach my $ac (split(/,/, $proteinForce)){
    $proteinForce{$ac} = 1;
  }
}

my %basicModifs;
foreach my $m (split(/,/, $basicModifList)){
  $basicModifs{$m} = 1;
}

my (%elemMass, %aaMass, %modifMassShift);
parseDefinitions($definitionFile);
my $cTermMass = $elemMass{H} + $elemMass{O};
my $nTermMass = $elemMass{H};
my $HpMass = $elemMass{'H+'};

# Parses files
my (%peptDescr, %prot, %spectra, $itemOrder, %GpeptDescr, %Gprot, %Gspectra);
our $file;
my $fileNum = 0;
foreach (@ARGV){
  foreach $file (glob($_)){
    print STDERR "Parsing $file\n" if ($verbose);
    if (!$parseAll && (($file =~ /-BLK\d+-/) || ($file =~ /-ACN\d+-/) || ($file =~ /-BSA-/))){
      print STDERR "Skipped\n" if ($verbose);
      next;
    }

    undef(%peptDescr);
    undef(%prot);
    undef(%spectra);

    if ($file =~ /\.gz$/){
      open(F, "gunzip -c $file |") || print STDERR "Warning, cannot open [$file]: $!";
    }
    else{
      open(F, $file) || print STDERR "Warning, cannot open [$file]: $!";
    }
    pidresParse(\*F, $fileNum++);
    close(F);

    unless (defined($notOnlyBestPept)){
      # Only keeps in %spectra the peptide giving the highest score
      foreach my $sp (keys(%spectra)){
	my @scores = sort {$a <=> $b} values(%{$spectra{$sp}{peptides}});
	my $maxScore = $scores[-1];
	foreach my $pm (keys(%{$spectra{$sp}{peptides}})){
	  if ($spectra{$sp}{peptides}{$pm} < $maxScore){
	    delete($spectra{$sp}{peptides}{$pm});
	    delete($peptDescr{$pm});

            # FIXME: inefficient deleting of match in protein
            foreach my $ac (keys(%prot)) {
              delete($prot{$ac}{matches}{$pm}) if defined ($prot{$ac}{matches}{$pm});
            }
	  }
	}
      }
    }

    # Selects and saves peptide/spectrum matches
    if ($noFileLevel){
      foreach my $sp (keys(%spectra)){
	$Gspectra{$sp} = $spectra{$sp};
      }
      foreach my $pm (keys(%peptDescr)){
	$GpeptDescr{$pm} = $peptDescr{$pm};
      }
      foreach my $ac (keys(%prot)){
	foreach my $match (keys(%{$prot{$ac}{matches}})){
	  $Gprot{$ac}{matches}{$match} = 1;
	}
      }
    }
    else{
      # Number of distinct peptides per protein that are above $minScore
      my (%numPept, %distinct, %protScore, %sphScore);
      foreach my $ac (keys(%prot)){
	my %pept;
	foreach my $match (keys(%{$prot{$ac}{matches}})){
	  my $zScore = $peptDescr{$match}{zScore};
	  my $pValue = $peptDescr{$match}{pValue};
	  my $peptide = $peptDescr{$match}{seq};
          my $specialModifOK = $peptDescr{$match}{specialModifOK};
	  if ($proteinForce{$ac} || ($specialModifOK && ($zScore >= $minScore) && ($pValue <= $maxPvalue) && (length($peptide) >= $minLen) && ($peptide =~ $qrValidAASeq))){
	    $pept{$peptide} = $zScore if (!defined($pept{$peptide}) || $zScore > $pept{$peptide});
	  }
	  if ($specialModifOK && ($zScore > $distinct{$ac}{$peptide})){
	    $distinct{$ac}{$peptide} = $zScore;
	  }
	}
	$numPept{$ac} = scalar(keys(%pept));
	foreach my $pept (keys(%pept)){
	  $sphScore{$ac} += $pept{$pept};
	}
	foreach my $score (values(%{$distinct{$ac}})){
	  $protScore{$ac} += $score;
	}
	#print STDERR "$ac has $numPept{$ac} distinct peptides, prot score $protScore{$ac} (".join('+',sort {$a<=>$b} values(%{$distinct{$ac}})).", sph score $sphScore{$ac} (".join('+',sort {$a<=>$b} values(%{$pept{$ac}})).")\n";
      }

      # Puts in global structure
      foreach my $ac (keys(%prot)){
	my $seqcov = 0;
	if ($saveScore < $defaultSaveScore){
	  # Single peptide hits are accepted, collect seq cov
	  my ($seq) = undef; ## GETSEQ($ac);
	  if ($seq){
	    my @seq = split(//, $seq);
	    foreach my $pept (keys(%{$distinct{$ac}})){
	      my $pos = -1;
	      while (($pos = index($seq, $pept, $pos+1)) != -1){
		for (my $i = $pos; $i < $pos+length($pept); $i++){
		  $seq[$i] = '*';
		}
	      }
	    }
	    my $count;
	    foreach my $char (@seq){
	      $count++ if ($char eq '*');
	    }
	    $seqcov = $count/length($seq);
	  }
	}
	if ($proteinForce{$ac} || (($protScore{$ac} >= $minProtScore) && (($numPept{$ac} >= $minNumPept) || (($sphScore{$ac} >= $saveScore) && ($seqcov >= $minSeqCovSPH))))){
	  #print STDERR "SPH $ac $file\n" if ($numPept{$ac} == 1);
	  foreach my $match (keys(%{$prot{$ac}{matches}})){
	    my $zScore = $peptDescr{$match}{zScore};
	    my $pValue = $peptDescr{$match}{pValue};
	    my $peptide = $peptDescr{$match}{seq};
            my $specialModifOK = $peptDescr{$match}{specialModifOK};
            my $containsSpecialModif = $peptDescr{$match}{containsSpecialModif};
	    if ($proteinForce{$ac} || ($specialModifOK && (!$containsSpecialModif || ($zScore >= $specialModifScore)) && ($zScore >= $outputScore) && ($pValue <= $maxPvalue) && (length($peptide) >= $minLen) && ($peptide =~ $qrValidAASeq))){
	      $Gprot{$ac}{matches}{$match} = 1;
	      $GpeptDescr{$match} = $peptDescr{$match};
	      $Gspectra{$peptDescr{$match}{cmpd}} = $spectra{$peptDescr{$match}{cmpd}};
	    }
	  }
	}
      }
    }
  }
}

# Global validation

# Detects proteins identified by the same set of peptides exactly
my (%pattern, %acToPattern);
foreach my $ac (keys(%Gprot)){
  my $pattern = join('|', sort({$a cmp $b} keys(%{$Gprot{$ac}{matches}})));
  $acToPattern{$ac} = $pattern;
  if ($pattern{$pattern}){
    # Another protein has the same set of queries
    #print STDERR "Protein $ac shares its peptides with $pattern{$pattern} [$pattern]\n";
    undef($Gprot{$ac}) if (defined($distinctProt));
  }
  else{
    $pattern{$pattern} = $ac;
  }
}

# Number of distinct peptides per protein that are above $minScore
my (%numPept, %distinct, %seqcov, %sphScore, %protScore);
foreach my $ac (keys(%Gprot)){
  my %pept;
  foreach my $match (keys(%{$Gprot{$ac}{matches}})){
    my $zScore = $GpeptDescr{$match}{zScore};
    my $pValue = $GpeptDescr{$match}{pValue};
    my $peptide = $GpeptDescr{$match}{seq};
    my $specialModifOK = $GpeptDescr{$match}{specialModifOK};
    if ($proteinForce{$ac} || ($specialModifOK && ($zScore >= $minScore) && ($pValue <= $maxPvalue) && (length($peptide) >= $minLen) && ($peptide =~ $qrValidAASeq))){
      $pept{$peptide} = $zScore if (!defined($pept{$peptide}) || $zScore > $pept{$peptide});
    }
    if ($specialModifOK && ($zScore > $distinct{$ac}{$peptide})){
      $distinct{$ac}{$peptide} = $zScore;
    }
  }
  $numPept{$ac} = scalar(keys(%pept));
  foreach my $pept (keys(%pept)){
    $sphScore{$ac} += $pept{$pept};
  }
  foreach my $score (values(%{$distinct{$ac}})){
    $protScore{$ac} += $score;
  }
  #print STDERR "$ac has $numPept{$ac} distinct peptides, prot score $protScore{$ac} (".join('+',sort {$a<=>$b} values(%{$distinct{$ac}})).", sph score $sphScore{$ac} (".join('+',sort {$a<=>$b} values(%{$pept{$ac}})).")\n";

  if ($saveScore < $defaultSaveScore){
    # Recomputes seqcov since now all samples are pooled and it may change. Saves in a hash to use for both big reds and final validation
    $seqcov{$ac} = 0;
    if ($saveScore < $defaultSaveScore){
      # Single peptide hits are accepted, collect seq cov
      my ($seq) = undef; ##GETSEQ($ac);
      if ($seq){
	my @seq = split(//, $seq);
	foreach my $pept (keys(%{$distinct{$ac}})){
	  my $pos = -1;
	  while (($pos = index($seq, $pept, $pos+1)) != -1){
	    for (my $i = $pos; $i < $pos+length($pept); $i++){
	      $seq[$i] = '*';
	    }
	  }
	}
	my $count;
	foreach my $char (@seq){
	  $count++ if ($char eq '*');
	}
	$seqcov{$ac} = $count/length($seq);
      }
    }
  }
}

# big reds computation, proteins identified based on the same peptides are considered equivalent
my %bigRedProt;
if ($minBigRed > 0){
  my %already;
  foreach my $ac (sort {-$protScore{$a} <=> -$protScore{$b}} keys(%Gprot)){
    if (($protScore{$ac} >= $minProtScore) && (($numPept{$ac} >= $minNumPept) || (($protScore{$ac} >= $saveScore) && ($seqcov{$ac} >= $minSeqCovSPH)))){
      my $numBigRed;
      foreach my $match (keys(%{$Gprot{$ac}{matches}})){
	my $peptide = $GpeptDescr{$match}{seq};
	if (!$already{$peptide}){
	  # First time we assign this peptide to a protein
	  $numBigRed++;
	  $already{$peptide} = $ac;
	}
	elsif ($acToPattern{$ac} eq $acToPattern{$already{$peptide}}){
	  # Already assigned but to a protein that is identified based on the same peptides exactly
	  $numBigRed++;
	}
      }
      $bigRedProt{$ac} = 1 if ($numBigRed >= $minBigRed);
    }
  }
}

# Starts printing
my $cmdLine = "pidresParser.pl".(defined($verbose)?' -verbose':'')." --savescore=$saveScore --outputscore=$outputScore --minscore=$minScore --maxpvalue=$maxPvalue --minnumpept=$minNumPept --minprotscore=$minProtScore --instrument=$instrument --minlen=$minLen";

my @time = localtime();
my $date = sprintf("%d-%02d-%02d", 1900+$time[5], 1+$time[4], $time[3]);
my $time = sprintf("%02d:%02d:%02d", $time[2], $time[1], $time[0]);
print <<end_of_xml;
<?xml version="1.0" encoding="ISO-8859-1"?>
  <idi:ProtSpectraIdentifications  xmlns:idi="namespace/ProtSpectra.html">
    <idi:OneSample>
      <idi:header>
        <idi:instrument>$instrument</idi:instrument>
        <idi:spectrumType>msms</idi:spectrumType>
        <idi:date>$date</idi:date>
        <idi:time>$time</idi:time>
        <idi:searchEngines>
          <idi:oneEngine>
            <idi:engineName>Phenyx</idi:engineName>
            <idi:engineVersion>$phenyxVersion</idi:engineVersion>
            <idi:searchDB dbName="$database" dbRelease="$dbRelease" dbReleaseDate="$dbReleaseDate"/>
            <idi:parentTol value="$searchParentTol" unit="$searchParentTolUnit"/>
            <idi:fragmentTol value="$searchFragmentTol" unit="$searchFragmentTolUnit"/>
            <idi:massType>$massType</idi:massType>
            <idi:cleavageEnzyme name="$enzymeName" missedCleavage="$nmc" cleavMode="$cleaveMode"/>
            <idi:fixedModifs>$fixedModifs</idi:fixedModifs>
            <idi:variableModifs>$variableModifs</idi:variableModifs>
            <idi:validationParams>
              <idi:autoExtraction><![CDATA[$cmdLine]]></idi:autoExtraction>
              <idi:saveScore>$saveScore</idi:saveScore>
              <idi:defaultSaveScore>$defaultSaveScore</idi:defaultSaveScore>
              <idi:minSeqCovSPH>$minSeqCovSPH</idi:minSeqCovSPH>
              <idi:outputScore>$outputScore</idi:outputScore>
              <idi:minScore>$minScore</idi:minScore>
              <idi:maximumPValue>$maxPvalue</idi:maximumPValue>
              <idi:minNumPept>$minNumPept</idi:minNumPept>
              <idi:minProtScore>$minProtScore</idi:minProtScore>
              <idi:minLen>$minLen</idi:minLen>
            </idi:validationParams>
          </idi:oneEngine>
        </idi:searchEngines>
$itemOrder      </idi:header>
    <idi:Identifications>
end_of_xml

# Selects and print peptide/spectrum matches
foreach my $ac (sort {-$protScore{$a} <=> -$protScore{$b}} keys(%Gprot)){
  if ($proteinForce{$ac} || (($protScore{$ac} >= $minProtScore) && (($numPept{$ac} >= $minNumPept) || (($sphScore{$ac} >= $saveScore) && ($seqcov{$ac} >= $minSeqCovSPH))) && (($minBigRed <= 0) || $bigRedProt{$ac}))){
    #print STDERR "SPH $ac $file\n" if ($numPept{$ac} == 1);
    print STDERR "$ac -------------------------------\n" if ($verbose);
    print "      <idi:OneProtein>\n      <idi:proteinId>$ac</idi:proteinId>\n      <idi:protSumScore>$protScore{$ac}</idi:protSumScore>\n";
    foreach my $match (sort {-$GpeptDescr{$a}{zScore} <=> -$GpeptDescr{$b}{zScore}} keys(%{$Gprot{$ac}{matches}})){
      my $zScore = $GpeptDescr{$match}{zScore};
      my $peptide = $GpeptDescr{$match}{seq};
      my $pValue = $GpeptDescr{$match}{pValue};
      my $specialModifOK = $GpeptDescr{$match}{specialModifOK};
      my $containsSpecialModif = $GpeptDescr{$match}{containsSpecialModif};
      if ($proteinForce{$ac} || ($specialModifOK && (!$containsSpecialModif || ($zScore >= $specialModifScore)) && (length($peptide) >= $minLen) && ($peptide =~ $qrValidAASeq) && ($zScore >= $outputScore) && ($pValue <= $maxPvalue))){
	my $modif = $GpeptDescr{$match}{modif};
	my $start = $GpeptDescr{$match}{start};
	my @modif = split(/:/, $modif);
	print STDERR "$peptide ($match)\n" if ($verbose);
        my $theoMass = getPeptideMass(pept=>$peptide, modif=>\@modif);
	my ($charge, $moz2) = getCorrectCharge($theoMass, $Gspectra{$GpeptDescr{$match}{cmpd}}{expMoz});
	print <<end_of_xml;
      <idi:OneIdentification>
        <idi:answer>
          <idi:sequence>$peptide</idi:sequence>
          <idi:modif>$modif</idi:modif>
          <idi:theoMass>$theoMass</idi:theoMass>
          <idi:charge>$charge</idi:charge>
          <idi:startPos>$start</idi:startPos>
          <idi:retentionTime>$Gspectra{$GpeptDescr{$match}{cmpd}}{rt}</idi:retentionTime>
        </idi:answer>
        <idi:source>
          <idi:file>$GpeptDescr{$match}{file}</idi:file>
          <idi:peptScore engine="Phenyx" pValue="$pValue">$zScore</idi:peptScore>
        </idi:source>
        <ple:peptide spectrumKey="$Gspectra{$GpeptDescr{$match}{cmpd}}{title}" xmlns:ple="namespace/PeakListExport.html">
$Gspectra{$GpeptDescr{$match}{cmpd}}{masslist}      </idi:OneIdentification>
end_of_xml
      }
    }
    print "      </idi:OneProtein>\n";
  }
}

print "    </idi:Identifications>
  </idi:OneSample>
</idi:ProtSpectraIdentifications>\n";

exit(0);


# --------------------------- subroutines -----------------------------------


sub pidresParse
{
  my ($F, $fileNum) = @_;
  my ($peptideKey, $sequence, $modif, $zScore, $zValue, $pValue, $charge, $ac, $match, $fm, $vm, $theoMass);

  my $shortFileName = basename($file);
  $shortFileName = (split(/\.pidres/, $shortFileName))[0];
  while (<$F>){

    if (/<dbmgr:origFastaFile>.*?([^\/]+)_v(\d+\.\d+)_(\d+)\..*<\/dbmgr:origFastaFile>/){
      # old format
      if (defined($database) && ($database ne $1)){
	die("Changed search database [$database, $1, $file]\n");
      }
      $database = $1;
      if (defined($dbRelease) && ($dbRelease ne $2)){
	die("Changed search database release [$dbRelease, $2, $file]\n");
      }
      $dbRelease = $2;
      if (defined($dbReleaseDate) && ($dbReleaseDate ne $3)){
	die("Changed search database release date [$dbReleaseDate, $3, $file]\n");
      }
      $dbReleaseDate = $3;
    }
    elsif (/<dbmgr:origFastaFile>(.+)<\/dbmgr:origFastaFile>/){
      # new format
      $database = $1;
    }
    elsif (/<dbmgr:release>(.+)<\/dbmgr:release>/){
      # new format
      $dbRelease = $1;
      $dbReleaseDate = $1.'01';
      $dbReleaseDate =~ s/\.//go;
    }
    elsif (/<peptTolerance.*unit="([^"]+)".*value="([^"]+)"/){
      if (defined($searchParentTolUnit) && ($searchParentTolUnit ne $1)){
	die("Changed search parent tolerance unit [$searchParentTolUnit, $1, $file]\n");
      }
      $searchParentTolUnit = $1;
      $searchParentTol = $2;
    }
    elsif (/<olavWideBaf:fragTolerance.*value="([^"]+)".*unit="([^"]+)"/){
      if (defined($searchFragmentTolUnit) && ($searchFragmentTolUnit ne $2)){
	die("Changed search parent tolerance unit [$searchFragmentTolUnit, $2, $file]\n");
      }
      $searchFragmentTolUnit = $2;
      $searchFragmentTol = $1;
    }
    elsif (/<Modifications>/){
      $fm = $vm = '';
    }
    elsif (/<OneModif.*mode="([^"]+)".*name="([^"]+)"/){
      if ($1 eq 'fixed'){
	$fm .= ',' if (length($fm) > 0);
	$fm .= $2;
      }
      else{
	$vm .= ',' if (length($vm) > 0);
	$vm .= $2;
      }
    }
    elsif (/<\/Modifications>/){
      if (defined($fixedModifs) && ($fixedModifs ne $fm)){
	die("Changed search fixed modifs [$fixedModifs; $fm; $file]\n");
      }
      $fixedModifs = $fm;
      if (defined($variableModifs) && ($variableModifs ne $vm)){
	print STDERR "Changed search variable modifs [$variableModifs; $vm; $file]\n" if (defined($verbose));
      }
      $variableModifs = $vm;
    }
    elsif (/<cleavageEnzyme.*name="([^"]+)".*missedCleavage="(\d+)".*cleavMode="([^"]+)"/){
      if (defined($enzymeName) && ($enzymeName ne $1)){
	die("Changed search enzyme [$enzymeName; $1; $file]\n");
      }
      $enzymeName = $1;
      $nmc = $2;
      $cleaveMode = $3;
    }

    # Peptides
    elsif (/<PeptideMatchDef .*key="([^"]+)">/){
      $peptideKey = $1;
    }
    elsif (/<sequence>(.+)<\/sequence>/){
      $sequence = $1;
    }
    elsif (/<modifAt>(.+)<\/modifAt>/){
      $modif = $1;
    }
    elsif (/<peptZScore>(.+)<\/peptZScore>/){
      $zScore = $1;
    }
    elsif (/<peptZvalue>(.+)<\/peptZvalue>/){
      $zValue = $1;
    }
    elsif (/<charge>(.+)<\/charge>/){
      $charge = $1;
    }
    elsif (/<\/PeptideMatchDef>/){
      my $pKey = "$fileNum-$peptideKey";
      $peptDescr{$pKey}{seq} = $sequence;
      $peptDescr{$pKey}{modif} = $modif;
      $peptDescr{$pKey}{zScore} = $zScore;
      $peptDescr{$pKey}{zValue} = $zValue;
      $peptDescr{$pKey}{charge} = $charge;
      my $specialModifOK = 1;
      my $containsSpecialModif = 0;
      foreach my $sm (split(/:/,$modif)){
        if (length($sm) > 1){
          if (!$basicModifs{$sm}){
            $containsSpecialModif = 1;
            if ($zScore < $specialModifScore){
              $specialModifOK = 0;
              last;
            }
          }
        }
      }
      $peptDescr{$pKey}{specialModifOK} = $specialModifOK;
      $peptDescr{$pKey}{containsSpecialModif} = $containsSpecialModif;
      my $cmpd = (split(/cmpd_/, $peptideKey))[1];
      $peptDescr{$pKey}{cmpd} = "$fileNum-$cmpd";
      $spectra{"$fileNum-$cmpd"}{peptides}{$pKey} = $zScore;
      $peptDescr{$pKey}{file} = $shortFileName;
    }

    # Proteins
    elsif (/<AC>(.*)<\/AC>/){
      $ac = $1;
    }
    elsif (/<PeptideMatch .*ref="([^"]+)">/){
      $match = "$fileNum-$1";
      $prot{$ac}{matches}{$match} = 1;
    }
    elsif (/<start>(.+)<\/start>/){
      $peptDescr{$match}{start} = $1+1; # We want first position is 1
    }
    elsif (/<pValue>(.+)<\/pValue>/){
      $peptDescr{$match}{pValue} = $1;
    }

    # Spectra
    elsif (/<ple:ItemOrder/){
      $itemOrder = $_;
      while (<$F>){
	$itemOrder .= $_;
	last if (/<\/ple:ItemOrder>/);
      }
    }
    elsif (/<ple:peptide/){
      my $compoundNumber;
      $internSequence++;
      if (/<ple:peptide .*compoundNumber="([^"]+)"/){
	$compoundNumber = $1;
      }
      else{
	$compoundNumber = $internSequence;
      }
      my $cmpdKey = "$fileNum-$compoundNumber";
      my $title;
      while (<$F>){
	if (/<ple:acquTime>(.+)<\/ple:acquTime>/){
	  $spectra{$cmpdKey}{rt} = $1;
	}
	else{
	  if (/<ple:PeptideDescr><!\[CDATA\[(.+)\]\]><\/ple:PeptideDescr>/){
	    $title = $1;
	    $spectra{$cmpdKey}{masslist} .= $_;
	  }
	  elsif (/<ple:ParentMass><!\[CDATA\[([.0-9]+) ([.0-9]+).*\]\]><\/ple:ParentMass>/){
	    $spectra{$cmpdKey}{expMoz} = $1;
	    $spectra{$cmpdKey}{intensity} = $2;
	    $spectra{$cmpdKey}{masslist} .= $_;
	  }
	  else{
	    $spectra{$cmpdKey}{masslist} .= $_ unless(defined($lightXML) && /^\d+\.\d+\s/);
	  }
	}
	last if (/<\/ple:peptide>/);
      }
      # To be compatible with Mascot, which loses spectrum orders we have to repeat the title and cannot use compoundNumber
      $spectra{$cmpdKey}{key} = $title; #($compoundNumber+1)."|$shortFileName";
      $spectra{$cmpdKey}{title} = $title;
   }
  }

} # pidresParse


sub dbConnect
{
  my ($db) = @_;

  return DBI->connect("dbi:Pg:dbname=$db", "biodbprod", "#4biodbprod", {PrintError=>0, RaiseError=>1, AutoCommit=>0}) 
    or die "Could not connnect to $db";

} # dbConnect


sub getPeptideMass
{
  my (%h) = @_;
  my ($pept, $modif) = ($h{pept}, $h{modif});
  die("No peptide given in getPeptideMass") unless (defined($pept));

  my $mass = 0.0;
  unless ($pept =~ /qrValidAASeq/){
    if (defined($modif)){
      # Applies all the mass shifts
      foreach (@$modif){
	$mass += $modifMassShift{$_} if (length($_) > 0);
      }
    }
    foreach (split(//, $pept)){
      $mass += $aaMass{$_};
    }
    $mass += $nTermMass+$cTermMass;
    return $mass;
  }
  else{
    # No defined mass
    print STDERR "no mass for peptide [$pept]\n";
    return -1.0 ;
  }

} # getPeptideMass


sub getCorrectCharge
{
  my ($mTheo, $mExp, $delta, $maxCharge) = @_;

  # Mass tolerance in Da
  $delta = $delta || 5;

  # Maximum charge
  $maxCharge = $maxCharge || 5;

  for (my $z = 1; $z <= $maxCharge; $z++){
    my $theoMoz = ($mTheo+$z*$HpMass)/$z;
    return ($z, $theoMoz) if (abs($theoMoz-$mExp) < $delta);
  }
  return (0,0);

} # getCorrectCharge


sub parseDefinitions
{
  my $fname = shift;

  use XML::Parser;
  open(F,$fname) || die("Cannot open definition file [$fname]\n");
  my $parser = new XML::Parser(Style => 'Stream');
  $parser->parse(\*F);
  close(F);

} # parseDefinitions


my ($modresName, $curChar);
sub StartTag
{
  my ($p, $el) = @_;

  if ($el eq 'oneElement'){
    $elemMass{$_{symbol}} = $_{monoisotopic};
  }
  elsif ($el eq 'oneAminoAcid'){
    $aaMass{$_{code1}} = $_{monoisotopic};
  }
  elsif ($el eq 'oneModRes'){
    $modresName = $_{name};
  }
  elsif ($el eq 'delta'){
    $modifMassShift{$modresName} = $_{monoisotopic};
    undef($modresName);
  }

  undef($curChar);

} # StartTag


sub Text
{
  $curChar .= $_;

} # Text


sub EndTag
{
  my($p, $el)= @_;

} # EndTag

sub get_val {
  my ($sect,$key) = @_;
  if ($config->exists(uc($sect),$key)) {
    return($config->val(uc($sect),$key));
  } else {
    return($config->val("DEFAULT",$key));
  }
}

sub get_params {
  return("PARAMETER SET $defaultSet:\n".
  sprintf("          minScore = %10s\n",$minScore).
  sprintf("       outputScore = %10s\n",$outputScore).
  sprintf("         maxPvalue = %10s\n",$maxPvalue).
  sprintf("  defaultSaveScore = %10s\n",$defaultSaveScore).
  sprintf("         saveScore = %10s\n",$saveScore).
  sprintf(" specialModifScore = %10s\n",$specialModifScore).
  sprintf("    basicModifList = %10s\n",$basicModifList).
  sprintf("      minSeqCovSPH = %10s\n",$minSeqCovSPH).
  sprintf("      minProtScore = %10s\n",$minProtScore).
  sprintf("        minNumPept = %10s\n",$minNumPept).
  sprintf("            minLen = %10s\n",$minLen).
  sprintf("         minBigRed = %10s\n",$minBigRed));
}
