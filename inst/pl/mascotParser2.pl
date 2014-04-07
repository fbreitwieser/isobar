#!/usr/bin/perl

# Mass spectrometry Perl program for extracting correct peptide matches from Mascot .dat files

# Copyright (C) 2007-2013 Jacques Colinge

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
#  Lazarettgasse 14
#  A-1090 Vienna, Austria
#  www.cemm.at

=head1 NAME

mascotParser2.pl - Extraction of reliable protein/peptide/spectrum matches from Mascot .dat files.

=head1 SYNOPSIS

mascotParser2.pl --defaultset=SETNAME [options] .dat files.

=head1 OPTIONS

Use mascotParser2.pl -h

=head1 DESCRIPTION

The script parses one or several Mascot .dat files to extract reliable peptide/spectrum matches and
outputs them in the .protSpectra.xml format. The .dat file(s) can be compressed (gzipped) files.

The selection of the peptide assignments is performed based on several thresholds applied to identifications
found in the .dat file(s):

=over 4

=item -help -h

=item -verbose

=item -parseall

Parses all the files, otherwise BSA, BLK, and ACN runs are skipped.

=item -notonlybestpept

Not only keeps the best peptide for each spectrum but all of them.

=item -distinctprot

Only outputs proteins with at least one specific peptide.

=item -nofilelevel

Global instead of file-based check of the selection criteria. Useful for considering several chromatographic/gel fractions as a single sample.

=item --proteinforce=AC-list

Let specific proteins pass through any threshold (but the basicscore) for QC or increased sensitivity, AC-list must be comma-separated. One use is to follow specific peptides in AP tags.

=item --minbigred=int

Number of specific peptides to be a big red (see Mascot documentation); only proteins having a big red peptide are selected.

=item --minscore=float

Minimum ion score (Mascot peptide score).

=item --minprotscore=float

Minimum protein score.

=item --minnumpept=int

Minimum number of distinct peptides per protein.

=item --savescore=float

Minimum peptide save ion score to bypass the minimum number of peptides required to be above minscore (allows single peptide hits with a higher threshold). Default is +infinity.

=item --minseqcovsph=float

Minimum protein sequence coverage for single peptide hits. Can be used to avoid gigantic proteins identified on the basis of a single peptide hit.

=item --minlen=int

Minimum peptide sequence length.

=item --basicscore=float

Minimum ion score to read a peptide from the .dat file (simple pre-filtering).

=item --outputscore=float

Minimum ion score to putput peptides for a selected protein.

=item --defaultset=['cid','cid_sensitive','hcd','qtof','qtof_sensitive']

Set the parameters to defaults according to a set of predifined values. The default values set not using this parameter are more conservative. The command line changes overwrite the defaults. Hence it is possible to take a default set and to change a specific parameter through the command line.

=item -lightXML

Flag to avoid outputing mass lists (much smaller XML). Default: yes.

=item -no-lightXML

Export mass lists.

=item -modifconv-file=string

File for conversion of modification strings to InsilicoSpectro equivalent. By default, the file modifconv.csv in the script directory is used.

=item -no-modifconv

Do not convert modifications.

=back

To be selected a peptide must have an ion score larger than minscore. The protein is finally selected if it has minnumpept peptides above minscore at least and a protein score better than minprotscore.

Alternatively, it is possible to specify savescore and then the condition on the number of peptides with score better than minscore is substitued by the following one: sum of the peptide scores better than minscore must be above savescore and a minimum sequence coverage of minseqcovsph must be reached. Sequence coverage is determined using all the peptides > outputscore. This alternative validation is mostly used for single peptide hits (SPHs).

In every case, as soon as a protein is selected, all the peptides above outputscore (usually < minscore) are reported.

During the parsing of the file, each spectrum is associated with the peptide that gives the best match, i.e. all multiple interpretations of a spectrum are lost in favor of the best one, unless -notonlybestpept flag is set. Moreover, all peptides with score less than the basic score (typically 5) are not read.

=head1 EXAMPLE

./mascotParser2.pl --defaultset=cid example.dat > test.protSpectra.xml

=head1 AUTHOR

Jacques Colinge

=cut

use strict;
use warnings;
no warnings 'uninitialized';
use FindBin qw/$Bin/;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Statistics::Lite qw(:all);
use Config::IniFiles;

my ($help, $verbose);

my $instrument = 'n/a';
my $modifconv_file = $Bin."/modifconv.csv";
my $configFile = $Bin."/mascotParser.ini";
my ($notOnlyBestPept, $distinctProt, $noFileLevel, $no_modifconv,
    $parseAll, $newDefaults, $hcdDefaults, $lightXML, $proteinForce, $defaultSet);

$lightXML=1;

my $config = Config::IniFiles->new(-file=>$configFile,-default=>'DEFAULT') 
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
  print STDERR "Use one of the available parameter sets using --defaultset=SETNAME:\n".
      "\t",join("\n\t",$config->Sections())."\n";
  exit(1);
}

my $basicScore       = $config->val(uc($defaultSet),'basicScore');
my $outputScore      = $config->val(uc($defaultSet),'outputScore');
my $minScore         = $config->val(uc($defaultSet),'minScore');
my $defaultSaveScore = $config->val(uc($defaultSet),'defaultSaveScore');
my $saveScore        = $config->val(uc($defaultSet),'saveScore');
my $minSeqCovSPH     = $config->val(uc($defaultSet),'minSeqCovSPH');
my $minProtScore     = $config->val(uc($defaultSet),'minProtScore');
my $minNumPept       = $config->val(uc($defaultSet),'minNumPept');
my $minLen           = $config->val(uc($defaultSet),'minLen');
my $minBigRed        = $config->val(uc($defaultSet),'minBigRed');


my $result = 
  GetOptions('help' => \$help,
             'h' => \$help,
             'parseall' => \$parseAll,
             'notonlybestpept' => \$notOnlyBestPept,
             'distinctprot' => \$distinctProt,
             'nofilelevel' => \$noFileLevel,
             'lightXML!' => \$lightXML,
             'proteinforce=s' => \$proteinForce,
             'minbigred=i' => \$minBigRed,
             'basicscore=f' => \$basicScore,
             'savescore=f' => \$saveScore,
             'defaultsavescore=f' => \$defaultSaveScore,
             'minseqcovsph=f' => \$minSeqCovSPH,
             'outputscore=f' => \$outputScore,
             'minscore=f' => \$minScore,
             'minnumpept=i' => \$minNumPept,
             'minprotscore=f' => \$minProtScore,
             'minlen=i' => \$minLen,
             'instrument=s' => \$instrument,
             'no-modifconv' => \$no_modifconv,
             'modifconv-file=s' => \$modifconv_file,
             'verbose' => \$verbose);

pod2usage(-verbose=>2, -exitval=>2) if (!$result);
pod2usage(get_params()) if (defined($help) && defined($defaultSet));
pod2usage(-verbose=>2, -exitval=>2) if (defined($help));
pod2usage(".dat files must be provided on command line!") if (!@ARGV || scalar(@ARGV) == 0);

if (defined($outputScore)) {
  pod2usage("outputscore (for hits of validated proteins) may not be higher than minscore (for validation)!") if ($outputScore > $minScore);
  pod2usage("outputscore (for hits of validated proteins) must be higher than basicscore (for pre-filtering)!") if ($outputScore < $basicScore);
} else {
  pod2usage("basicscore (for pre-filtering) must be higher than minscore (for validation)!") if ($basicScore > $minScore);
}


$outputScore = $minScore if (!defined($outputScore));

my %proteinForce;
if (defined($proteinForce)){
  foreach my $ac (split(/,/, $proteinForce)){
    $proteinForce{$ac} = 1;
  }
}

# Charge format conversion
my %charge = (
	      '1+,2+,and3+' => '1,2,3',
	      '1+,2+and3+' => '1,2,3',
	      '1+' => '1',
	      '2+' => '2',
	      '2+and3+' => '2,3',
	      '2+,and3+' => '2,3',
	      '3+' => '3',
	      '4+' => '4'
	     );

# Modifications conversion (Macot mod_file to InSilicoSpectro insilicodef.xml names)
my $modifConv = read_csvhash($modifconv_file);

# Parses files
my (%cmpd, %prot, %query, %Gcmpd, %Gprot);
my ($database, $comment, $dbRelease, $mascotVersion, $dbReleaseDate);
my ($searchFragmentTol, $searchFragmentTolUnit, $searchParentTol, $searchParentTolUnit);
my ($fixedModifs, $variableModifs, $nmc, $enzymeName, $cleaveMode);
my $massType = 'monoisotopic';
my (%elemMass, %aaMass, %modifMassShift, @nonTermFixedModif, @fixedModifLocation, @fixedModif, @variableModif);
my ($nTermFixedModif, $cTermFixedModif, $cTermMass, $nTermMass, $HpMass);
my $qrValidAASeq = qr/^[ACDEFGHIJKLMNOPQRSTUVWY]+$/;
our $file;
my $fileNum = 0;
foreach (@ARGV){
  foreach $file (glob($_)){
    print STDERR "Parsing $file\n" if ($verbose);
    if (!$parseAll && (($file =~ /-(BLK\d*|ACN\d+|BSA|CLN)-/))){
      print STDERR "Skipped\n" if ($verbose);
      next;
    }

    undef(%cmpd);
    undef(%prot);
    undef(%query);
    undef(@fixedModif);
    undef(@variableModif);
    undef(@nonTermFixedModif);
    undef($cTermFixedModif);
    undef($nTermFixedModif);

    if ($file =~ /\.gz$/){
      open(F, "gunzip -c $file |") || print STDERR "Warning, cannot open [$file]: $!";
    }
    else{
      open(F, $file) or die "Warning, cannot open [$file]: $!";
    }
    mascotParse(\*F, $fileNum++);
    close(F);

    unless (defined($notOnlyBestPept)){
      # Only keeps in %prot and %query the best score for each peptide
      foreach my $query (keys(%query)){
	my @scores = sort {$a <=> $b} values(%{$query{$query}{ac}});
	my $bestScore = $scores[-1];
	foreach my $ac (keys(%{$query{$query}{ac}})){
	  if ($prot{$ac}{queries}{$query}{score} < $bestScore){
	    undef($prot{$ac}{queries}{$query});
	    undef($query{$query}{ac}{$ac});
	  }
	}
      }
    }

    # Selects and saves peptide/spectrum matches, i.e. no file level validation
    if ($noFileLevel){
      foreach my $ac (keys(%prot)){
	foreach my $query (keys(%{$prot{$ac}{queries}})){
	  $Gprot{$ac}{queries}{$query} = $prot{$ac}{queries}{$query};
	  $Gcmpd{$query} = $cmpd{$query};
	}
      }
    }
    else{
      # File level validation

      # Number of distinct peptides per protein that are above $minScore
      my (%numPept, %distinct, %sphScore, %protScore);
      foreach my $ac (keys(%prot)){
	my %pept;
	foreach my $query (keys(%{$prot{$ac}{queries}})){
	  my $score = $prot{$ac}{queries}{$query}{score};
	  my $peptide = $prot{$ac}{queries}{$query}{pept};
	  if ($proteinForce{$ac} || (($score >= $minScore) && (length($peptide) >= $minLen) && ($peptide =~ $qrValidAASeq))){
	    $pept{$peptide} = $score if (!defined($pept{$peptide}) || $score > $pept{$peptide});
	  }
	  if ($score > $distinct{$ac}{$peptide}){
	    $distinct{$ac}{$peptide} = $score;
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

      # Puts in global structure if good enough
      foreach my $ac (keys(%prot)){
	my $seqcov = 0;
	if ($saveScore < $defaultSaveScore){
	  # Single peptide hits are accepted, collect seq cov
	  my ($seq) = undef; # disfunct: Need database to get sequence
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
	  foreach my $query (keys(%{$prot{$ac}{queries}})){
	    my $peptide = $prot{$ac}{queries}{$query}{pept};
	    my $score = $prot{$ac}{queries}{$query}{score};
	    if ($proteinForce{$ac} || ((length($peptide) >= $minLen) && ($peptide =~ $qrValidAASeq) && ($score >= $outputScore))){
	      $Gprot{$ac}{queries}{$query} = $prot{$ac}{queries}{$query};
	      $Gcmpd{$query} = $cmpd{$query};
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
  ## FB: Perl warns that sort is not numeric and does not sort, thus
  ##     I removed sorting for now as it does not change output. Pattern: O-12345
  ## my $pattern = join('|', sort({$a <=> $b} keys(%{$Gprot{$ac}{queries}})));
  my $pattern = join('|', keys(%{$Gprot{$ac}{queries}}));
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
  foreach my $query (keys(%{$Gprot{$ac}{queries}})){
    my $score = $Gprot{$ac}{queries}{$query}{score};
    my $peptide = $Gprot{$ac}{queries}{$query}{pept};
    if ($proteinForce{$ac} || (($score > $minScore) && (length($peptide) >= $minLen) && ($peptide =~ $qrValidAASeq))){
      $pept{$peptide} = $score if (!defined($pept{$peptide}) || $score > $pept{$peptide});
    }
    if ($score > $distinct{$ac}{$peptide}){
      $distinct{$ac}{$peptide} = $score;
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
      my ($seq) = undef; # disfunct: needs database connection to get sequence
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
    if (($protScore{$ac} >= $minProtScore) && (($numPept{$ac} >= $minNumPept) || (($sphScore{$ac} >= $saveScore) && ($seqcov{$ac} >= $minSeqCovSPH)))){
      my $numBigRed;
      foreach my $query (keys(%{$Gprot{$ac}{queries}})){
	my $peptide = $Gprot{$ac}{queries}{$query}{pept};
	if (!$already{$peptide}){
	  # First time we assign this peptide to a protein
	  $numBigRed++;
	  $already{$peptide} = $ac;
	}
	elsif ($acToPattern{$ac} eq $acToPattern{$already{$peptide}}){
	  # Already assigne but to a protein that is identified based on the same peptides exactly
	  $numBigRed++;
	}
      }
      $bigRedProt{$ac} = 1 if ($numBigRed >= $minBigRed);
    }
  }
}

# Starts printing
my $cmdLine = "mascotParser.pl".(defined($verbose)?' -verbose':'')." --basicscore=$basicScore --savescore=$saveScore --minseqcovsph=$minSeqCovSPH --outputscore=$outputScore --minscore=$minScore --minnumpept=$minNumPept --minprotscore=$minProtScore --instrument=$instrument --minlen=$minLen";

my @time = localtime();
my $date = sprintf("%d-%02d-%02d", 1900+$time[5], 1+$time[4], $time[3]);
my $time = sprintf("%02d:%02d:%02d", $time[2], $time[1], $time[0]);
my $lxml = defined($lightXML) ? 'yes' : 'no';
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
            <idi:engineName>Mascot</idi:engineName>
            <idi:engineVersion>$mascotVersion</idi:engineVersion>
            <idi:searchDB dbName="$database" dbRelease="$dbRelease" dbReleaseDate="$dbReleaseDate"/>
            <idi:parentTol value="$searchParentTol" unit="$searchParentTolUnit"/>
            <idi:fragmentTol value="$searchFragmentTol" unit="$searchFragmentTolUnit"/>
            <idi:massType>$massType</idi:massType>
            <idi:cleavageEnzyme name="$enzymeName" missedCleavage="$nmc" cleavMode="$cleaveMode"/>
            <idi:fixedModifs>$fixedModifs</idi:fixedModifs>
            <idi:variableModifs>$variableModifs</idi:variableModifs>
            <idi:validationParams>
              <idi:autoExtraction><![CDATA[$cmdLine]]></idi:autoExtraction>
              <idi:basicScore>$basicScore</idi:basicScore>
              <idi:saveScore>$saveScore</idi:saveScore>
              <idi:defaultSaveScore>$defaultSaveScore</idi:defaultSaveScore>
              <idi:minSeqCovSPH>$minSeqCovSPH</idi:minSeqCovSPH>
              <idi:outputScore>$outputScore</idi:outputScore>
              <idi:minScore>$minScore</idi:minScore>
              <idi:minNumPept>$minNumPept</idi:minNumPept>
              <idi:minProtScore>$minProtScore</idi:minProtScore>
              <idi:minLen>$minLen</idi:minLen>
              <idi:lightXML>$lxml</idi:lightXML>
            </idi:validationParams>
          </idi:oneEngine>
        </idi:searchEngines>
	<ple:ItemOrder xmlns:ple="namespace/PeakListExport.html">
	  <ple:item type="mass"/>
	  <ple:item type="intensity"/>
	  <ple:item type="charge"/>
	</ple:ItemOrder>
      </idi:header>
    <idi:Identifications>
end_of_xml

# Selects and print peptide/spectrum matches
foreach my $ac (sort {-$protScore{$a} <=> -$protScore{$b}} keys(%Gprot)){
  if ($proteinForce{$ac} || (($protScore{$ac} >= $minProtScore) && (($numPept{$ac} >= $minNumPept) || (($sphScore{$ac} >= $saveScore) && ($seqcov{$ac} >= $minSeqCovSPH))) && (($minBigRed <= 0) || $bigRedProt{$ac}))){
    #print STDERR "global SPH $ac $file\n" if ($numPept{$ac} == 1);
    print STDERR "$ac -------------------------------\n" if ($verbose);
    print "      <idi:OneProtein>\n      <idi:proteinId>$ac</idi:proteinId>\n      <idi:protSumScore>$protScore{$ac}</idi:protSumScore>\n";

    ## FB: Perl warns that sort is not numeric and does not sort, thus
    ##     I removed the sorting for now as it does not change output. Pattern: O-12345
    ## foreach my $query (sort {$a <=> $b} keys(%{$Gprot{$ac}{queries}}))
    foreach my $query (keys(%{$Gprot{$ac}{queries}})){
      my $peptQry = $Gprot{$ac}{queries}{$query};
      my $peptide = $peptQry->{pept};
      my $peptide_length_ge_min = (length($peptide) >= $minLen);
      my $peptide_isvalid = ($peptide =~ $qrValidAASeq);
      my $peptide_score_ge_outputScore =  ($peptQry->{score} >= $outputScore);
      if ($proteinForce{$ac} || ($peptide_length_ge_min && $peptide_isvalid && $peptide_score_ge_outputScore)){
	my ($modifStr, @modif);
	convertMascotModif($peptide, $peptQry->{modif}, \$modifStr, \@modif);
	if ($modifStr){
	  my $theoMass = getPeptideMass(pept=>$peptide, modif=>\@modif);
	  my ($charge, $moz2) = getCorrectCharge($theoMass, $Gcmpd{$query}{expMoz});
	  my $pValue = sprintf "%.2e", exp(-log(10)*$peptQry->{score}/10);

	  print STDERR "$peptide ($query)\n" if ($verbose);
	  print <<end_of_xml;
      <idi:OneIdentification>
        <idi:answer>
          <idi:sequence>$peptide</idi:sequence>
          <idi:modif>$modifStr</idi:modif>
          <idi:theoMass>$theoMass</idi:theoMass>
          <idi:charge>$charge</idi:charge>
          <idi:startPos>$peptQry->{start}</idi:startPos>
          <idi:aaBefore>$peptQry->{aaBefore}</idi:aaBefore>
          <idi:aaAfter>$peptQry->{aaAfter}</idi:aaAfter>
          <idi:retentionTime>$Gcmpd{$query}{rt}</idi:retentionTime>
        </idi:answer>
        <idi:source>
          <idi:file>$Gcmpd{$query}{file}</idi:file>
          <idi:peptScore engine="Mascot" pValue="$pValue">$Gprot{$ac}{queries}{$query}{score}</idi:peptScore>
        </idi:source>
        <ple:peptide spectrumKey="$Gcmpd{$query}{title}" xmlns:ple="namespace/PeakListExport.html">
        <ple:PeptideDescr><![CDATA[$Gcmpd{$query}{title}]]></ple:PeptideDescr>
        <ple:ParentMass><![CDATA[$Gcmpd{$query}{expMoz} $Gcmpd{$query}{intensity} $Gcmpd{$query}{charge}]]></ple:ParentMass>
        <ple:peaks><![CDATA[
end_of_xml
	  if (defined($lightXML)){
	    print "]]></ple:peaks>
        </ple:peptide>
      </idi:OneIdentification>\n";
	  }
	  else{
	    print "$Gcmpd{$query}{massList}]]></ple:peaks>
        </ple:peptide>
      </idi:OneIdentification>\n";
	  }
	}
	else{
	  print STDERR "Cannot use peptide [$peptide] because one modification has no InSilicoSpectro equivalent [$Gprot{$ac}{queries}{$query}{modif}]\n" if ($verbose);
	}
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


sub convertMascotModif
{
  my ($pept, $mascotModif, $modifStr, $vModif) = @_;

  # First set the variable modification as found by Mascot
  my @mascotModif = split(//, $mascotModif);
  for (my $i = 0; $i < @mascotModif; $i++){
    $mascotModif[$i] = $variableModif[$mascotModif[$i]];
  }

  # To compute masses, we only need the variable modifications as Mascot includes the fixed modifications in the AA masses or at the n- c-term masses
  @$vModif = @mascotModif;

  # Second, we localize the fixed terminal mofifications
  if (defined($nTermFixedModif)){
    $mascotModif[0] = $fixedModif[$nTermFixedModif];
  }
  if (defined($cTermFixedModif)){
    $mascotModif[-1] = $fixedModif[$cTermFixedModif];
  }

  # Third, we complement with the fixed AA modifs to get the modif string
  my @pept = split(//,$pept);
  for (my $i = 0; $i < @pept; $i++){
    foreach my $m (@nonTermFixedModif){
      if ($pept[$i] =~ $fixedModifLocation[$m]){
	$mascotModif[$i+1] = $fixedModif[$m];
	last;
      }
    }
  }
  $$modifStr = join(':',@mascotModif);

} # convertMascotModif


sub mascotParse
{
  my ($F, $fileNum) = @_;

  my $shortFileName = basename($file);
  my $qryMatch;
  my %affectedAcs;

  $shortFileName = (split(/\.dat/, $shortFileName))[0]; # Get rid of the Mascot result file extension if any
  while (<$F>){
    s/[\r\n]//go;

    if (/^Content-Type: +application\/x\-Mascot; +name="parameters"/){
      parseParameters($F);
    }
    if (/^Content-Type: +application\/x\-Mascot; +name="header"/){
      parseHeader($F);
    }
    if (/^Content-Type: +application\/x\-Mascot; +name="masses"/){
      parseMasses($F);
    }
    elsif (/^q(\d+)_p(\d+)=(.*)/){
      my $query = $1;
      $qryMatch = $2;
      my ($nmc, $mass, $delta, $nIons, $pept, $nUsed1, $modif, $score, $ionSeries, $nUsed2, $nUsed3, @part) = split(/[,;:](?![^,;:"]+")/, $3);
      if ($score >= $basicScore){
	for (my $i = 0; $i < @part/5; $i++){
	  my $ac = $part[5*$i];
	  $ac =~ s/"//go;
	  if ($ac eq "0"){
	    print STDERR "$i > $_\n";
	    exit(0);
	  }
	  my $qKey = "$fileNum-$query";
	  if (!defined($prot{$ac}{queries}{$qKey}{score}) || ($prot{$ac}{queries}{$qKey}{score} < $score)){ # test in case the same spectrum matches several peptides in the same DB entry
	    $prot{$ac}{queries}{$qKey} = {
	    	score => $score,
	        pept => $pept,
	        modif => $modif,
	        start => $part[5*$i+2],
            aaBefore => '?',
            aaAfter => '?' };
	    $query{$qKey}{ac}{$ac} = $score;
	    # remember the position of the protein
        $affectedAcs{$ac} = $i;
	  }
	}
      } else {
	  undef( %affectedAcs );
	  }
	}
	elsif (/^q(\d+)_p(\d+)_terms=(.*)/){
	  # flanking residues
	  my $query = $1;
	  if ( $qryMatch == $2 && %affectedAcs ) {
	    my $qKey = "$fileNum-$query";
	    my @terms = split(/[,:]/, $3);
	    while ( my ($ac, $index) = each %affectedAcs) {
          $prot{$ac}{queries}{$qKey}{aaBefore} = $terms[2*$index];
          $prot{$ac}{queries}{$qKey}{aaAfter} = $terms[2*$index+1];
        }
      }
	  undef( %affectedAcs );
	  undef( $qryMatch );
    }
    elsif (/^qmass(\d+)=(.+)/){
      $cmpd{"$fileNum-$1"}{expMass} = $2;
    }
    elsif (/^qexp(\d+)=(.+),(.+)/){
      $cmpd{"$fileNum-$1"}{expMoz} = $2;
    }
    elsif (/^qintensity(\d+)=(.+)/){
      $cmpd{"$fileNum-$1"}{intensity} = $2;
    }
    elsif (/^Content\-Type: +application\/x\-Mascot; +name="query(\d+)"/){
      parseOneExpSpectrum("$fileNum-$1", $shortFileName, $1, $F);
    }
  }

} # mascotParse


sub parseParameters
{
  my $F = shift;

  while (<$F>){
    last if (index($_, '--gc0p4Jq0M2Yt08jU534c0p') == 0);

    s/[\r\n]//go;
    if (/^DB=(.+)/){
      if (defined($database) && ($database ne $1)){
	die("Changed search database [$database, $1, $file]\n");
      }
      $database = $1;
    }
    elsif (/^COM=(.+)/){
      $comment = $1;
    }
    elsif (/^TOL=(.+)/){
      $searchParentTol = $1;
    }
    elsif (/^TOLU=(.+)/){
      if (defined($searchParentTolUnit) && ($searchParentTolUnit ne $1)){
	die("Changed search parent tolerance unit [$searchParentTolUnit, $1, $file]\n");
      }
      $searchParentTolUnit = $1;
    }
    elsif (/^ITOL=(.+)/){
      $searchFragmentTol = $1;
    }
    elsif (/^ITOLU=(.+)/){
      if (defined($searchFragmentTolUnit) && ($searchFragmentTolUnit ne $1)){
	die("Changed search fragment tolerance unit [$searchParentTolUnit, $1, $file]\n");
      }
      $searchFragmentTolUnit = $1;
    }
    elsif (/^CLE=(.+)/){
      if (defined($enzymeName) && ($enzymeName ne $1)){
	die("Changed search enzyme [$enzymeName, $1, $file]\n");
      }
      $enzymeName = $1;
      $cleaveMode = 'normal';
    }
    elsif (/^PFA=(.+)/){
      $nmc = $1
    }
  }

} # parseParameters


sub parseMasses
{
  my $F = shift;

  $variableModif[0] = '';
  my ($Hmass, $eMass);
  while (<$F>){
    last if (index($_, '--gc0p4Jq0M2Yt08jU534c0p') == 0);

    s/[\r\n]//go;
    if (/^([A-Z])=(.+)/){
      $aaMass{$1} = $2;
    }
    elsif (/^Hydrogen=(.+)/){
      $Hmass = $1;
    }
    elsif (/^Electron=(.+)/){
      $eMass = $1;
    }
    elsif (/^C_term=(.+)/){
      $cTermMass = $1;
    }
    elsif (/^N_term=(.+)/){
      $nTermMass = $1;
    }
    elsif (/^FixedMod([1-9])=([^,]+),(.+)/){
      if ($no_modifconv) {
        $fixedModif[$1-1] = replaceColon($3);
      } else {
        $fixedModif[$1-1] = $modifConv->{trim($3)};
        if (!$fixedModif[$1-1]){
          print STDERR "ERROR: No InSilicoSpectro equivalent defined for fixed modification [$3].\n";
          print STDERR "       in modifconv file (--modifconv-file=$modifconv_file)\n.";
          print STDERR "       Either defined the modification in the modifconv file, or set --no-modifconv\n";
          exit 1;
        }
      }
      $modifMassShift{$fixedModif[$1-1]} = $2;
    }
    elsif (/^FixedModResidues([1-9])=(.+)/){
      $fixedModifLocation[$1-1] = $2;
      if ($fixedModifLocation[$1-1] eq 'N_term'){
	$nTermFixedModif = $1-1;
      }
      elsif ($fixedModifLocation[$1-1] eq 'C_term'){
	$cTermFixedModif = $1-1;
      }
      else{
	push(@nonTermFixedModif, $1-1);
	my $pattern = $fixedModifLocation[$1-1];
	$fixedModifLocation[$1-1] = qr/[$pattern]/;
      }
    }
    elsif (/^delta([1-9])=([^,]+),(.+)/){
      if ($no_modifconv) {
        $fixedModif[$1-1] = replaceColon($3);
      } else {
        $variableModif[$1] = $modifConv->{trim($3)};
        if (!$variableModif[$1]){
          print STDERR "ERROR: No InSilicoSpectro equivalent defined for fixed modification [$3].\n";
          print STDERR "       in modifconv file (--modifconv-file=$modifconv_file)\n.";
          print STDERR "       Either defined the modification in the modifconv file, or set --no-modifconv\n";
          exit 1;
        }
      }
      $modifMassShift{$variableModif[$1]} = $2;
    }
  }
  $HpMass = $Hmass-$eMass;
  my $fm = join(',',@fixedModif);
  if (defined($fixedModifs) && ($fixedModifs ne $fm)){
    die("Changed search fixed modifs [$fixedModifs; $fm; $file]\n");
  }
  $fixedModifs = $fm;
  my $vm = join(',',@variableModif[1..$#variableModif]);
  if (defined($variableModifs) && ($variableModifs ne $vm)){
    die("Changed search variable modifs [$variableModifs; $vm; $file]\n");
  }
  $variableModifs = $vm;

} # parseMasses

sub replaceColon {
  $_[0] =~ s/:/_/g;
  return ($_[0]);
}

sub parseHeader
{
  my $F = shift;

  while (<$F>){
    last if (index($_, '--gc0p4Jq0M2Yt08jU534c0p') == 0);

    s/[\r\n]//go;
    if (/^version=(.+)/){
      if (defined($mascotVersion) && ($mascotVersion ne $1)){
        die("Changed Mascot version [$mascotVersion, $1, $file]\n");
      }
      $mascotVersion = $1;
    }
    elsif (/^release=.+_v(\d+\.\d+)_(\d+)/){
      if (defined($dbRelease) && ($dbRelease ne $1)){
        die("Changed search database release [$dbRelease, $1, $file]\n");
      }
      if (defined($dbReleaseDate) && ($dbReleaseDate ne $2)){
        die("Changed search database release date [$dbReleaseDate, $2, $file]\n");
      }
      $dbRelease = $1;
      $dbReleaseDate = $2;
    }
  }

} # parseHeader


sub parseOneExpSpectrum
{
  my ($queryNum, $shortFileName, $query, $F) = @_;

  my ($massList, $charge, $rt, $title);
  while (<$F>){
    last if (index($_, '--gc0p4Jq0M2Yt08jU534c0p') == 0);

    s/[\r\n]//go;
    if (/charge=(.+)/){
      $charge = $1;
      $charge =~ s/\s//go;
      $charge = $charge{$charge};
      $charge = '2,3' if (!$charge);
    }
    elsif (/title=(.+)/){
      $title = $1;
      $title =~ s/%2d/-/go;
      $title =~ s/%2e/./go;
    }
    elsif (/rtinseconds=(.+)/){
      $rt = $1;
      if ($rt =~ /([\d\.]+)\-([\d\.]+)/){
        # The retention time is given as a range, take the average
        $rt = ($1+$2)*0.5;
      }
    }
    elsif (/Ions1=(.+)/){
      my @part = split(/,/, $1);
      undef($massList);
      foreach my $peak (@part){
        my ($mass, $intensity) = split(/:/, $peak);
        $massList .= "$mass $intensity ?\n";
      }
    }
  }

  $cmpd{$queryNum}{charge} = $charge;
  $cmpd{$queryNum}{rt} = $rt;
  $cmpd{$queryNum}{title} = $title;
  $cmpd{$queryNum}{key} = $title; # Mascot loses the spectrum order and $query is the order after sorting "$query|$shortFileName";
  $cmpd{$queryNum}{massList} = $massList;
  $cmpd{$queryNum}{file} = $shortFileName;

} # parseOneExpSpectrum


sub dbConnect
{
  my ($db) = @_;

  return DBI->connect("dbi:Pg:dbname=$db", "biodbprod", "#4biodbprod", {PrintError=>0, RaiseError=>1, AutoCommit=>0});

} # dbConnect


sub getPeptideMass
{
  my (%h) = @_;
  my ($pept, $modif) = ($h{pept}, $h{modif});
  die("No peptide given in getPeptideMass") unless (defined($pept));

  my $mass = 0.0;
  unless ($pept =~ /InSilicoSpectro::InSilico::AASequence::qrValidAASeq/){
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

sub read_csvhash {
   my ($csvhash_file) = @_;
   my %csvHash;
   print STDERR "csvHash file $csvhash_file:\n" if $verbose;
   open (CSVHASH, $csvhash_file) or die "Could not open $csvhash_file:  $!";
   while (<CSVHASH>) {
      if (!/^\s*#/ && !/^\s*$/) {
         chomp();
         my ($tmpvar1, $tmpvar2) = split(/\t+/, $_);
         if (defined($tmpvar1) && defined($tmpvar2)) {
            $csvHash{trim($tmpvar1)} = trim($tmpvar2);
            print STDERR "[".trim($tmpvar1)."] => [".trim($tmpvar2)."]\n" if $verbose;
         }
      }
   }
   close(CSVHASH);
   return(\%csvHash);
} #read_csvhash

# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
  my $string = shift;
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
}

sub get_params {
  return(
"Validation Principle: \n".
"            1.) Proteins must be validated \n".
"                 - with [->minNumPept,normally 2] peptide hits [-> minScore] or \n".
"                 - a high scoring single peptide hit [-> saveScore]\n".
"            2.) Peptides of validated proteins may be accepted with lower scores
                 [-> outputScore]".
"    
    PARAMETER SET $defaultSet:\n".
  sprintf("          minScore = %10s\tminimum score for validation (requires [minNumPept] hits with score >= minScore)\n",$minScore).
  sprintf("         saveScore = %10s\tminimum score for single peptide hits of proteins\n",$saveScore).
  sprintf("       outputScore = %10s\tminimum score for peptides of validated proteins\n\n",$outputScore).

  sprintf("        basicScore = %10s\tpre-filtering score -has no influence results, must be lower than the others\n",$basicScore).
  sprintf("  defaultSaveScore = %10s\tignore\n",$defaultSaveScore).
  sprintf("      minSeqCovSPH = %10s\tminimum sequence coverage for single peptide hits [not working]\n",$minSeqCovSPH).
  sprintf("      minProtScore = %10s\tminimum protein score\n",$minProtScore).
  sprintf("        minNumPept = %10s\tminimum number of peptides\n",$minNumPept).
  sprintf("            minLen = %10s\tminimum length of pepties\n",$minLen).
  sprintf("         minBigRed = %10s\tminimum number of 'big reds' per protein\n",$minBigRed));
}
