#!/usr/bin/perl
# Creation date : 2010-09-29
# Last modified : Fri 20 Jan 2012 10:26:15 AM CET

# Module        : tab2xls.pl
# Purpose       : converts csv files to XLS format
# Usage         : tab2xls.pl newfile.xls 'Sheet 1=file1.csv'
# Licence       : based on example files distributed with Spreadsheet::WriteExcel
#                 by John McNamara. reverse('©'), Dec 2010, Florian P Breitwieser
# Contact       : Florian P Breitwieser <fbreitwieser@cemm.oeaw.ac.at>

use strict;
use warnings;
use Spreadsheet::WriteExcel;
use File::Basename;

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

my $wbheader = $workbook->add_format();
$wbheader->set_color('white');
$wbheader->set_align('center');
#$wbheader->set_align('vcenter');
$wbheader->set_fg_color('green');

my $fmt_centeracross = $workbook->add_format(center_across=>1);
$fmt_centeracross->set_color('white');
$fmt_centeracross->set_pattern();
$fmt_centeracross->set_fg_color('green');

my $fmt_red = $workbook->add_format(color=>'red');
my $fmt_green= $workbook->add_format(color=>'blue');
my $fmt_gray= $workbook->add_format(color=>'gray');

my $wbcenter = $workbook->add_format();
$wbcenter->set_align('center');

for (my $i=0; $i <= $#ARGV; ++$i) {
    print STDERR " processing $ARGV[$i]\n";
    # Add a worksheet
    my ($name,$file) = split(/=/,$ARGV[$i]);
    
	my $props; 
	($name,$props) = get_props($name);

    if (!defined $file) {
    	$file = $name;
	    ($name) = getname(fileparse($name,".csv"));
    	#$name = "Sheet $i";
    }

    open(F,"<",$file) or die "Could not open $file: $!";
    
    my $header = <F>;
    my @header = split("\t",$header);
	 chomp(@header);
    my %header = map { $_ => 1 } @header;

    my $worksheet = $workbook->add_worksheet($name);
    $worksheet->add_write_handler(qr[\w], \&store_string_widths);

    my $freeze_col = (defined($props->{'freeze_col'})? $props->{'freeze_col'}:0);
    $worksheet->freeze_panes(1, $freeze_col); # 1 row
    my $col = 0;
    my $row = 0;

    for my $i (0..$#header) {
	    $header[$i] = trim($header[$i]);
        write_col($worksheet,$row,$col,$header[$i],1);
        ++$col;
    }
    $row = 1;
    
    while(my $line=<F>)    {
    	chomp($line);
	    my @data = split("\t",$line);
    	my $col = 0;
	    for (my $i=0; $i<scalar(@data); ++$i) {
    	    my $field = $data[$i];
            write_col($worksheet,$row,$col,trim($field),0);
	  	    ++$col;
    	}
    	++$row;
    }

	if (defined $props->{'autofilter'}) {
		$worksheet->autofilter(0,0,$row-1,$col-1);
	}

    close(F);
    
    # Run the autofit after you have finished writing strings to the workbook.
    autofit_columns($worksheet);
    
}

$workbook->close;


sub get_props {
	my ($string) = @_;
	my %def;
	if (my ($def,$f) = ($string =~ /^:(.*):(.*)$/)) {
    
    for my $s (split /,/,$def) {
      my ($key,$val) = split / /, $s;
      $val = 1 if !defined $val;
      print STDERR "key $key = $val\n";
      $def{$key} = $val;
    }
		$string = $f;
	}
	return ($string,\%def);
}


sub write_col {
    my ($worksheet,$row,$col,$field,$is_header) = @_;
    if (my ($def,$f) = ($field =~ /^:(.*):(.*)$/)) {
#		my %def = map {$_ => 1} split /,/,$def;
#		if (defined $def{'centeracross'}) {
		if ($def eq 'centeracross') {
		    $worksheet->write($row,$col,$f,$fmt_centeracross);
	    } else {
            die "no known col property";
        }
    } else {
        if ($is_header) {
    	    $worksheet->write($row,$col,getname($field),$wbheader);
        } else {
    	    if ($field eq 'TRUE') { $worksheet->write($row,$col,$field,$fmt_green);
            } elsif ($field eq '0') { $worksheet->write($row,$col,$field,$fmt_gray);
        	} else {
	       	    $worksheet->write($row,$col,$field);
    	    }
        }
    }
}

sub getname {
    my ($field) = @_;
    my @parts = map(ucfirst,split(/[\._]/,$field));
    return(join(" ",@parts));
}

sub trim {
	my ($string) = @_;
	$string =~ s/^"//;
	$string =~ s/"$//;
	return $string;
}

sub autofit_columns {

    my $worksheet = shift;
    my $col       = 0;

    for my $width (@{$worksheet->{__col_widths}}) {

        $worksheet->set_column($col, $col, $width) if $width;
        $col++;
    }
}

sub store_string_widths {

    my $worksheet = shift;
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

    # Ignore numbers
    return if $token =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/;

    # Ignore various internal and external hyperlinks. In a real scenario
    # you may wish to track the length of the optional strings used with
    # urls.
    return if $token =~ m{^[fh]tt?ps?://};
    return if $token =~ m{^mailto:};
    return if $token =~ m{^(?:in|ex)ternal:};


    # We store the string width as data in the Worksheet object. We use
    # a double underscore key name to avoid conflicts with future names.
    #
    my $old_width    = $worksheet->{__col_widths}->[$col];
    my $string_width = string_width($token);

    if (not defined $old_width or $string_width > $old_width) {
        # You may wish to set a minimum column width as follows.
        #return undef if $string_width < 10;

        $worksheet->{__col_widths}->[$col] = $string_width;
    }


    # Return control to write();
    return undef;
}

sub string_width {
    my $width = (0.9 * length $_[0]) + 2;
    if ($width>25) {$width=25;}
    return $width;
}
