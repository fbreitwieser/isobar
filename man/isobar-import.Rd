\name{isobar-import}

\alias{isobar-import}

% readIBSpectra
\alias{readIBSpectra}
\alias{readIBSpectra,character,character,character-method}
\alias{readIBSpectra,character,character,missing-method}
\alias{readIBSpectra,character,data.frame,missing-method}
\alias{readIBSpectra,character,data.frame,character-method}
\alias{readIBSpectra,character,character-method}
\alias{read.mzid}

\title{Loading data into IBSpectra objects using readIBSpectra}

\description{
        Read ibspectra-csv files and peaklist files as an
	IBSpectra object of type 'type' (see \code{\link{IBSpectra}},
	e.g. iTRAQ4plexSpectra or TMT6plexSpectra). If
	peaklist.file is missing, it is assumed that id.file contains
	intensity and m/z columns for the reporter tags.
}
      
\usage{
\S4method{readIBSpectra}{character,character}(type,id.file)

# reads id file
\S4method{readIBSpectra}{character,character,character}(
                 type, id.file, peaklist.file, sep = "\t", mapping.file
                 = NULL, mapping = c(quantification.spectrum = "hcd",
                 identification.spectrum = "cid"), id.file.domap =
                 NULL, identifications.format = NULL, decode.titles =
                 FALSE, ...)

# reads peaklist file
\S4method{readIBSpectra}{character,data.frame,character}(
                 type, id.file, peaklist.file, annotate.spectra.f =
                 NULL, peaklist.format = NULL, scan.lines = 0,
                 fragment.precision = NULL, fragment.outlier.prob =
                 NULL, ...)
}

\arguments{
  \item{type}{
    Name of class of new IBSpectra object: 
    \code{\link{iTRAQ4plexSpectra}}, \code{\link{iTRAQ8plexSpectra}},
    \code{\link{TMT2plexSpectra}}, \code{\link{TMT6plexSpectra}}, or \code{\link{TMT10plexSpectra}}
  }
  \item{id.file}{
    Database search results file in \code{ibspectra.csv} or
   \code{mzIdentML} format. See \code{identifications.format}. See the vignette for
   information on converting Mascot dat and Phenyx pidres files into
   ibspectra format.
  }
  \item{peaklist.file}{
    Peaklist file, typically in MGF format, see \code{peaklist.format}.
    MGF must be centroid!
  }
%  \item{proteinGroupTemplate}{
%    When having technical or biological repeats: First a template
%    protein group is created which uses information from all runs,
%    then this template is applied.
%    It should increase comparability across runs.
%  }
  \item{mapping.file}{
    If defined, spectum titles from the peaklist file are linked to
    the identifications via this file. This can be used when running HCD
    runs for quantification and CID runs for identification. See
    Koecher et al., 2009 for details.
  }
  \item{mapping}{
    Named character vector defining the names of columns
    in mapping.file. The names must be 'peaklist' and 'id', and
    the values must correspond to colnames of the mapping files.
  }
%  \item{mapping.file.readopts}{ 
%    Read options for \code{\link{read.table}} when reading files specified
%    in \code{mapping.file}.
%  }
  \item{id.file.domap}{
    When using HCD-CID or a method akin and every spectrum is used for 
    identification, the ID result files  of the HCD run can be
    specfied in \code{id.file.domap}. Then, the results are merged after
    mapping the identification results.
  }
  \item{annotate.spectra.f}{
    Function which changes or annotates the spectra feature data
    before it is written to IBpectra object. This can be used to calculate 
    and threshold additional scores, for example localization scores of post-
    translational modifications such as Delta Score (filterSpectraDeltaScore) or
    PhosphoRS site localization probabilities (annotateSpectraPhosphoRS).
  }
  \item{peaklist.format}{
    "mgf" (Mascot Generic format) or "mcn" (iTracker Machine Readable
    output). When NULL, it detects the format on file name extension.
  }
  \item{identifications.format}{
    "ibspectra.csv" or "mzid"  (PSI MzIdentML format). When NULL, file
    format is guessed based on extension.
  }
  \item{fragment.precision}{
    Fragment precision for extraction of reporter tags: for each tag
    and spectrum the m/z-intensity pair with it's mass closest to the
    known reporter tag mass is extracted within the window
    true_mass +/- fragment.precision/2.
  }
  \item{fragment.outlier.prob}{
    Fragment outlier probability filter: After all m/z-intensity pairs
    have been extracted, those pairs with the fragment.outlier.prob/2
    most unprecise m/z values are filtered out.
  }
  \item{decode.titles}{
    Boolean. Decode spectrum titles in identification file using \code{\link{URLdecode}}. 
    When extracting the  DAT file from Mascot web interface, the spectrum titles are
    encoded - \%20 instead of space, etc. Set decode.titles to TRUE to
    map these titles to the unescaped MGF titles.
  }
  \item{scan.lines}{
    Read files sequentially scan.lines lines at a time. Can help in
    case of memory issues, set to 10000 or higher, for example.
  }
  \item{sep}{sep argument of read.table}
  \item{\dots}{
    Further arguments handed down to \code{initialize}.
  }
}


\author{Florian P. Breitwieser, Jacques Colinge}

\seealso{
  \link{ProteinGroup},
  \link{IBSpectra},
  \link{isobar-preprocessing},
  \link{isobar-analysis},
  \link{isobar-plots}
}

\examples{
data(ibspiked_set1)

# get identifier for Ceruplasmin proteins
ceru.acs <- protein.g(proteinGroup(ibspiked_set1),"CERU")
# create a smaller ibspectra w/ only Ceruplasmins
ib.ceru <- subsetIBSpectra(ibspiked_set1,protein=ceru.acs,direction="include")

# write it to a file
tf <- tempfile("isobar")
write.table(as.data.frame(ib.ceru),sep="\t",file=tf,quote=FALSE)

# read it again into an IBSpectra object
ib.ceru2 <- readIBSpectra("iTRAQ4plexSpectra",tf,identifications.format="ibspectra")
ib.ceru2

unlink(tf)
}
