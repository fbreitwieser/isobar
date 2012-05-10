#!/usr/bin/perl
# Creation date : 2010-09-29
# Last modified : Wed 09 May 2012 05:59:39 PM CEST

# Module        : tab2xls.pl
# Purpose       : converts csv files to XLS format
# Usage         : tab2xls.pl newfile.xls 'Sheet 1=file1.csv'
# Licence       : based on example files distributed with Spreadsheet::WriteExcel
#                 by John McNamara. reverse('©'), Dec 2010, Florian P Breitwieser
# Contact       : Florian P Breitwieser <fbreitwieser@cemm.oeaw.ac.at>

use strict;
use warnings;
no warnings 'uninitialized';
use Spreadsheet::WriteExcel;
use File::Basename;
use Data::Dumper;

my $delim="\t";

# Check for valid number of arguments
if (($#ARGV < 0)) {
	print("Usage: tab2xls.pl newfile.xls ':PROPERTIES:Sheet 1=file1.csv'\n\n",
          "  Valid properties:\n",
          "    freeze_col N       freeze pane to column N\n",
          "    autofilter         apply autofilter on all columns\n"),
};

my $xls_file = $ARGV[0];
if (@ARGV == 1) {
  $xls_file =~ s/.csv/.xls/;
} else {
  shift;
}

print STDERR "writing XLS file $xls_file\n";

# Create a new Excel workbook
my $workbook = Spreadsheet::WriteExcel->new($xls_file);
$workbook->set_properties(
  comments=>'Created with isobar R package and Spreadsheet::WriteExcel');
#my $green = $workbook->set_custom_color(20, 44, 122, 45);
my $maroon = $workbook->set_custom_color(41, 158, 30, 30);
#my @colors = (
#  $workbook->set_custom_color(41, 158, 30, 30),
#  $workbook->set_custom_color(42, 60, 189, 60),
#  $workbook->set_custom_color(43, 60, 60, 189)
#);
my @colors = (
  $workbook->set_custom_color(42, 55, 126, 184),
  $workbook->set_custom_color(43, 77, 175, 74),
  $workbook->set_custom_color(41, 228, 26, 28),
  $workbook->set_custom_color(44, 152, 78, 163),
  $workbook->set_custom_color(45, 255, 127, 0),
  $workbook->set_custom_color(46, 255, 255, 51)
);


my %header_prop = (
  color => 'white',
  text_wrap => 1,
  align => 'center');

my $fmt_centeracross = $workbook->add_format(%header_prop,(center_across=>1,fg_color=>$colors[2]));

my %color_formats;
my $row_limit = 65536;
my @header_formats;

my $wbcenter = $workbook->add_format();
$wbcenter->set_align('center');

for (my $file_i=0; $file_i <= $#ARGV; ++$file_i) {
  print STDERR " processing $ARGV[$file_i]\n";
  # Add a worksheet

	my ($file,$props) = get_props($ARGV[$file_i],":");

  print STDERR "$file\n";

  my $name = $props->{'name'};
  if (!defined $name) { ($name) = getname(fileparse($file,".csv")); }

  open(F,"<",$file) or die "Could not open $file: $!";
    
  my $header = <F>;
  my @header = split("\t",$header);
  chomp(@header);
  my %header = map { $_ => 1 } @header;
  my ($worksheet,$row);
  my $n_worksheets = 0;
   
  while(my $line=<F>)    {
    if (!defined $row || $row == $row_limit) {
      $n_worksheets += 1;
      $worksheet = $workbook->add_worksheet($name.($n_worksheets > 1? " $n_worksheets" : ""));
      $worksheet->add_write_handler(qr[\w], \&store_string_widths);
      $worksheet->freeze_panes(1,(defined($props->{'freeze_col'})? $props->{'freeze_col'}:0)); # 1 row
      for my $i (0..$#header) {
        $header[$i] = trim($header[$i]);
        $header_formats[$i] = $workbook->add_format(%header_prop,fg_color=>$colors[$file_i]);
        write_col($worksheet,0,$i,getname($header[$i]),$header_formats[$i]);
      }
      $row = 1;
    }
    chomp($line);
    my @data = split("\t",$line);
    my $col = 0;
    my ($data0,$rowprops) = get_props($data[0],"#");
    $data[0] = $data0;
    if (defined $rowprops->{'color'}||defined $rowprops->{'bottomborder'}) {
      $worksheet->set_row($row,undef,colorfmt($rowprops->{'color'},$rowprops->{'bottomborder'}));
    }
    for (my $i=0; $i<scalar(@data); ++$i) {
      my $field = $data[$i];
      write_col($worksheet,$row,$col,trim($field));
      ++$col;
    }
    ++$row;
    if (defined $props->{'autofilter'} && $row == $row_limit) { $worksheet->autofilter(0,0,$row-1,$col-1); }
  }


  if (defined $props->{'autofilter'}) { $worksheet->autofilter(0,0,$row-1,$#header); }
  close(F);
    
  # Run the autofit after you have finished writing strings to the workbook.
  autofit_columns($worksheet);
    
}

$workbook->close;


sub get_props {
	my ($string,$sep) = @_;
  if (!defined $sep) { $sep = ":"; }
	my %def;
	if (my ($def,$f) = ($string =~ /^${sep}(.*)${sep}(.*)$/)) {
    
    for my $s (split /,/,$def) {
      my ($key,$val) = split /=/, $s;
      $val = 1 if !defined $val;
      $def{$key} = $val;
    }
		$string = $f;
	}
	return ($string,\%def);
}


sub write_col {
  my ($worksheet,$row,$col,$field,$format) = @_;

	my ($f,$props) = get_props($field,"@");
  #if (my ($def,$f) = ($field =~ /^|(^|*)|(.*)$/)) { ## only check for centeracross column property
#  if (my ($def,$f) = ($field =~ /^:([^:]+):(.*)$/)) {       ## this creates problem w/ modification column  (format: :::etc::)
  if ($props->{'centeracross'} eq 'centeracross') {
      $worksheet->write($row,$col,$f,$fmt_centeracross);
  } else {
    if ($field eq 'TRUE') { $format=colorfmt('green'); }
    elsif ($field eq '0') { $format=colorfmt('gray'); }
    $worksheet->write($row,$col,$f,$format);
  }

  if (defined $props->{'comment'}) {
    $worksheet->write_comment($row,$col,$props->{'comment'});
  }
}

sub colorfmt {
  my($color,$bottom) = @_;
  $bottom=0 unless defined $bottom;
  $color=0 unless defined $color;
  if (!defined $color_formats{$color}) {
    $color_formats{$color}{$bottom} = $workbook->add_format(color=>$color,bottom=>$bottom);
  }
  return($color_formats{$color}{$bottom});
}

sub getname {
    my ($field) = @_;
    $field =~ s/###/\n/;
    my @parts = map(ucfirst,split(/[\._]/,$field));
    return(join(" ",@parts));
}

sub trim {
	my ($string) = @_;
  return($string) unless defined $string;
	$string =~ s/^"//;
	$string =~ s/"$//;
	return $string;
}

sub autofit_columns {

    my $worksheet = shift;
    my $col       = 0;
    my $max_header = 4;

    for (my $i=0;$i < @{$worksheet->{__col_widths}}; ++$i) {
      my $width = ${$worksheet->{__col_widths}}[$i];
      my $header_width = ${$worksheet->{__header_widths}}[$i];
      if ($width && $header_width > $width) {
        $header_formats[$i]->set_rotation(90);
#        print STDERR "$i: $header_width\t$width\n";
        $max_header = $header_width if ($header_width > $max_header);
      }
      $worksheet->set_column($col, $col, $width) if $width;
      $col++;
    }
    $worksheet->set_row(0,$max_header*5);
}

sub store_string_widths {

    my $worksheet = shift;
    my $row       = $_[0];
    my $col       = $_[1];
    my $token     = $_[2];

    # Ignore some tokens that we aren't interested in.
    return if not defined $token;       # Ignore undefs.
    return if $token eq '';             # Ignore blank cells.
    return if ref $token eq 'ARRAY';    # Ignore array refs.
    return if $token =~ /^=/;           # Ignore formula
    return if $token =~ /^Isotope Im/;  # Ignore Isotope Imp
    return if $token =~ /^Analysis Properties/;  
    return if $token =~ /^Class Labels/;
#    return if $ro

    # Ignore numbers
    #return if $token =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/;
    my $is_number = $token =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/;

    # Ignore various internal and external hyperlinks. In a real scenario
    # you may wish to track the length of the optional strings used with
    # urls.
    return if $token =~ m{^[fh]tt?ps?://};
    return if $token =~ m{^mailto:};
    return if $token =~ m{^(?:in|ex)ternal:};


    # We store the string width as data in the Worksheet object. We use
    # a double underscore key name to avoid conflicts with future names.
    #
    my $string_width = string_width($token);
    if ($is_number && $string_width > 8) { $string_width = 8; }
    if ($row > 0) {
      my $old_width    = $worksheet->{__col_widths}->[$col];

      if (not defined $old_width or $string_width > $old_width) {
        # You may wish to set a minimum column width as follows.
        #return undef if $string_width < 10;

        $worksheet->{__col_widths}->[$col] = $string_width;
      }
    } else {
      # define the header width
      $worksheet->{__header_widths}->[$col] = $string_width;
    }


    # Return control to write();
    return undef;
}

sub string_width {
    my $max_width = 0;
    foreach (split(/\n/,$_[0])) {
      my $width = (0.9 * length $_) + 2;
      $max_width = $width if ($width>$max_width);
    }
    if ($max_width>25) {$max_width=25;}
    return $max_width;
}
