#!/usr/bin/perl
# Creation date : 2010-09-29
# Last modified : Thu 15 Aug 2013 11:41:48 PM CEST

# Module        : tab2xls.pl
# Purpose       : converts csv files to XLS format
# Usage         : tab2xls.pl newfile.xls 'Sheet 1=file1.csv'
# Licence       : based on example files distributed with Spreadsheet::WriteExcel
#                 by John McNamara. reverse('©'), Dec 2010, Florian P Breitwieser
# Contact       : Florian P Breitwieser <fbreitwieser@cemm.oeaw.ac.at>

use strict;
use warnings;
no warnings 'uninitialized';
use Excel::Writer::XLSX;
use File::Basename;

my $delim="\t";
my %conditional_formats;

# Check for valid number of arguments
if (($#ARGV < 0)) {
	print("Usage: tab2xls.pl newfile.xls ':PROPERTIES:Sheet 1=file1.csv'\n\n",
          "  Valid properties:\n",
          "    freeze_col N       freeze pane to column N\n",
          "    autofilter         apply autofilter on all columns\n"),
};

my ($xls_file,$props) = get_props($ARGV[0],":");
if (@ARGV == 1) {
  $xls_file =~ s/.csv/.xlsx/;
} else {
  shift;
}

print STDERR "writing XLS file $xls_file\n";

# Create a new Excel workbook
my $workbook = Excel::Writer::XLSX->new($xls_file);
$workbook->set_properties(
  title => $props->{'title'},
  author => $props->{'author'},
  subject => $props->{'subject'},
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
  text_wrap => 1,
  align => 'center');

my $fmt_centeracross = $workbook->add_format(align=>'center',valign=>'vcenter',center_across=>1,underline=>1);

my %color_formats;
my $row_limit = 65536;
my @header_formats;

my $wbcenter = $workbook->add_format();
$wbcenter->set_align('center');

for (my $file_i=0; $file_i <= $#ARGV; ++$file_i) {
  print STDERR " processing $ARGV[$file_i]\n";
  # Add a worksheet

	my ($file,$props) = get_props($ARGV[$file_i],":");

  #print STDERR "$file\n";

  my $name = $props->{'name'};
  #if (!defined $name) { ($name) = getname(fileparse($file,".csv")); }
  if (!defined $name) { ($name) = fileparse($file,".csv"); }

  open(F,"<",$file) or die "Could not open $file: $!";
    
  my $header = <F>;
  my @header = split("\t",$header);
  chomp(@header);
  my %header = map { $_ => 1 } @header;
  %conditional_formats = ();
  my ($worksheet,$row);
  my $n_worksheets = 0;
   
  while(my $line=<F>)    {
    if (!defined $row || $row == $row_limit) {
      $n_worksheets += 1;
      $worksheet = $workbook->add_worksheet(substr($name,0,25).($n_worksheets > 1? " $n_worksheets" : ""));
      $worksheet->add_write_handler(qr[\w], \&store_string_widths);
      my $fg_color = ($file_i <= $#colors)? "white" : "black";
      my $bg_color = ($file_i <= $#colors)? $colors[$file_i] : "white";
      $worksheet->set_tab_color($bg_color);
      if (defined $props->{'header'}) { $worksheet->set_header($props->{'header'}); }
      if (defined $props->{'footer'}) { $worksheet->set_footer($props->{'footer'}); }
      $worksheet->freeze_panes(1,(defined($props->{'freeze_col'})? $props->{'freeze_col'}:0)); # 1 row

      for my $i (0..$#header) {
        $header_formats[$i] = $workbook->add_format(%header_prop,
          color=>$fg_color,fg_color=>$bg_color);
      }
      write_header($worksheet,0,\@header_formats,@header);
      $row = 1;
    }
    #chop($line);
    my @data = split("\t",$line);
    my $col = 0;
    my ($data0,$rowprops) = get_props($data[0],"#");
    $data[0] = $data0;
    #print STDERR "data /vs/ header elems in line $row: $#data /vs/ $#header\n";
    while (scalar @data < scalar @header) {
      #print STDERR "push in next line\n";
      my $line2 = <F>; #chomp($line2);
      my @data2 = split("\t",$line2);
      $data[$#data] = $data[$#data].(shift @data2);
      @data = (@data,@data2);
    }
    if (scalar @data != scalar @header) {
      print STDERR "Bad CSV in row $row: $#data /vs/ $#header:\n";
      print "###### HEADER: ######\n",join("\n",@header),"\n\n";
      print "###### LINE $row #####\n",join("\n",@data),"\n";
      die();
    }
    chomp($data[$#data]);

    if (defined $rowprops->{'color'}||defined $rowprops->{'bottomborder'}) {
      $worksheet->set_row($row,undef,colorfmt($rowprops->{'color'},$rowprops->{'bottomborder'}));
    }
    if (defined $rowprops->{'level'}) {
      ## TODO: LEVEL
      #$worksheet->set_row($row,undef,undef
    }

    write_row($worksheet,$row,undef,$fmt_centeracross,\@data);

    ++$row;
    if (defined $props->{'autofilter'} && $row == $row_limit) { $worksheet->autofilter(0,0,$row-1,$col-1); }
  }


  if (defined $props->{'autofilter'}) { $worksheet->autofilter(0,0,$row-1,$#header); }
  close(F);
    
  # Run the autofit after you have finished writing strings to the workbook.
  autofit_columns($worksheet);
  set_conditional_formatting($worksheet,$row,\%conditional_formats);
    
}

$workbook->close;

sub set_conditional_formatting {
  my ($worksheet,$nrow,$cond_formats) = @_;

  while (my ($col,$format) = each %$cond_formats) {
    $worksheet->conditional_formatting(1,$col,$nrow,$col,
      {type=>$format,
      mid_color=>'#FFFFFF'})
  }
}

sub write_header {
  my ($worksheet,$row,$format_ref,@data) = @_;
  my ($merge_from,$do_merge,$merge_val,$mmerge_val);
  for (my $i=0; $i<scalar(@data); ++$i) {
    ($do_merge,$merge_val) = write_col($worksheet,$row,$i,trim($data[$i]),$format_ref->[$i]);
    if ($do_merge) {
      $merge_from = $i unless defined($merge_from);
      $mmerge_val = $merge_val unless defined $mmerge_val;
    } elsif (defined($merge_from)) {
      $worksheet->merge_range($row,$merge_from,$row,$i-1,
                              $mmerge_val,$format_ref->[$i]);
      $merge_from = undef;
      $mmerge_val = undef;
    }
  }
  if (defined($merge_from)) {
    $worksheet->merge_range($row,$merge_from,$row,$#data,
      $mmerge_val,$format_ref->[0]);
  }
}


sub write_row {
  my ($worksheet,$row,$format,$format_merge,$data) = @_;
  my ($merge_from,$do_merge,$merge_val,$mmerge_val);
  for (my $i=0; $i<scalar(@$data); ++$i) {
    ($do_merge,$merge_val) = write_col($worksheet,$row,$i,trim($data->[$i]),$format);
    if ($do_merge) {
      $merge_from = $i unless defined($merge_from);
      $mmerge_val = $merge_val unless defined $mmerge_val;
    } elsif (defined($merge_from)) {
      $worksheet->merge_range($row,$merge_from,$row,$i-1,
                              $mmerge_val,$format_merge);
      $merge_from = undef;
      $mmerge_val = undef;
    }
  }
  if (defined($merge_from)) {
    $worksheet->merge_range($row,$merge_from,$row,scalar(@$data)-1,
      $mmerge_val,$format_merge);
  }
}

sub get_props {
	my ($string,$sep) = @_;
  if (!defined $sep) { $sep = ":"; }
	my %def;
	if (my ($def,$f) = ($string =~ /^${sep}(.*)${sep}(.*)$/s)) {
    
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

  my $props;
  ($field,$props) = get_props($field,"@");
  if (defined $props->{'centeracross'}) {
    return (1,$field);
  } 
  if ($field eq 'TRUE') { $format=colorfmt('green'); }
  elsif ($field eq 'FALSE') { $format=colorfmt('lightgray'); }
  elsif ($field eq '0') { $format=colorfmt('gray'); }

  $field = getname($field) if ($row == 0);

  if (defined $props->{'link'}) {
    $worksheet->write_url($row,$col,$props->{'link'},$field,$format);
  } else {
    $worksheet->write($row,$col,$field,$format);
  }

  if (defined $props->{'comment'} && !($props->{'comment'} =~ /^[ \n]*$/)) {
    $props->{'comment'} =~ s/NULL//g;
    $props->{'comment'} =~ s/\n\s*\n+/\n/g;
    return if $props->{'comment'} =~ /^ *$/;
    #print STDERR "comment: [".$props->{'comment'}."]\n";
    $worksheet->write_comment($row,$col,$props->{'comment'},visible=>0);
  }

  if (defined $props->{'conditional_formatting'}) {
    if ($row == 0) {
      $conditional_formats{$col} = $props->{'conditional_formatting'};
    } else {
      $worksheet->conditional_formatting($row,$col,{type=>$props->{'conditional_formatting'}});
    }
  }
  return(0);
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
    return unless defined $worksheet->{__col_widths}; 

    for (my $i=0;$i < @{$worksheet->{__col_widths}}; ++$i) {
      my $width = ${$worksheet->{__col_widths}}[$i];
      my $header_width = ${$worksheet->{__header_widths}}[$i];
      if ($width && $header_width > $width) {
        $header_formats[$i]->set_rotation(90);
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
