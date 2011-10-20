\name{isobar-import}

\alias{isobar-import}

% readIBSpectra
\alias{readIBSpectra}
\alias{readIBSpectra,character,character,character-method}
\alias{readIBSpectra,character,character,missing-method}
\alias{readIBSpectra,character,character-method}
\alias{read.mzid}

\title{Loading data into IBSpectra objects using readIBSpectra}

\description{
        Read ibspectra-csv files and peaklist files as an
      IBSpectra object of type 'type' (see \code{\link{IBSpectra}},
      e.g. iTRAQ4plexSpectra or TMT6plexSpectra). If
      peaklist.file is missing, it assumes that id.file contains
      intensity and m/z information for the reporter tags.

}
      
\usage{
\S4method{readIBSpectra}{character,character}(type,id.file)
\S4method{readIBSpectra}{character,character,character}(type,id.file,
peaklist.file,proteinGroupTemplate=NULL,mapping.file=NULL,
mapping=c(peaklist=1,id=2),mapping.file.readopts = list(header = TRUE,
stringsAsFactors = FALSE, sep = ","),peaklist.format=NULL,id.format=NULL,
fragment.precision=NULL,fragment.outlier.prob=NULL,revert.titles=FALSE)
}



\arguments{
  \item{type}{
    \code{\link{iTRAQ4plexSpectra}}, \code{\link{iTRAQ8plexSpectra}},
    \code{\link{TMT2plexSpectra}}, or \code{\link{TMT6plexSpectra}}
  }
  \item{id.file}{
Database search results file in \code{ibspectra.csv} or
\code{mzIdentML} format. See \code{id.format}. See the vignette for
  information on converting Mascot dat and Phenyx pidres files into
  ibspectra format.
  }
  \item{peaklist.file}{
    Peaklist file, see \code{peaklist.format}.
  }
  \item{proteinGroupTemplate}{
When having technical or
biological repeats: First a template protein group is created
which uses information from all runs, then this template is
applied. It should increase comparability across runs.
}
  \item{mapping.file}{
If defined, spectum titles from the peaklist file are linked to
the identifications via this file. This can be used when running short HCD runs for
quantification and long CID runs for identification. See
Koecher et al., 2009 for details. }
  \item{mapping}{
mapping columns of mapping.file.
}
  \item{mapping.file.readopts}{
read options for read.table
  }
  \item{peaklist.format}{
"mgf" (Mascot Generic format) or "mcn" (iTracker Machine Readable
 output). When NULL, it detects the format on file name extension.
}
  \item{id.format}{
"ibspectra.csv" or "mzid"  (PSI MzIdentML format). When NULL, file
  format is guessed based on extension.
}

\item{fragment.precision}{
    Fragment precision for import of isobaric masses.
}
\item{fragment.outlier.prob}{
    Fragment outlier probability filter.
}

\item{revert.titles}{
    Unescape spectrum titles in identification file. Use with Mascot which might transform
    Spaces into \%20, etc.
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
\dontrun{
data(ibspiked_set1)
write.table(as.data.frame(ibspiked_set1),sep="\t",file="ib1.tmp.ibspectra.csv")
ib1 <- readIBSpectra("iTRAQ4plexSpectra","ib1.tmp.ibspectra.csv")
ib1
ibspiked_set1
file.remove("ib1.tmp.ibspectra.csv")
}
}
