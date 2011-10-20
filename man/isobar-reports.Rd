\name{isobar-reports}

\alias{isobar-reports}
\alias{initialize.report}

\title{Isobar reports}

\description{
  Generation of LaTeX and XLS reports is helped with functions
  which facilitate the gathering of relevant information and creation
  of tikz plots. \code{initialize.report} is the main function establishing 
  an environment with all neccessary data and computations.
}

\usage{
  initialize.report(properties.file = "properties.conf",
                    args = commandArgs(TRUE), env = .GlobalEnv)
}

\arguments{
  \item{properties.file}{File which holds the parameters for data analysis and report generation.
    It is written in standard R code and defines peaklists, identification files, significance levels,
    etc. See \code{system.file("report","properties.conf",package="isobar")} for 
    a full set of possible values. It is parsed using \code{\link{sys.source}}.}
  \item{args}{Additional parameters which overrids those in \code{properties.file}. }
  \item{env}{Environment to which the results are written.}
}

\details{
  The directory \code{inst} in the isobar installation directory 
  \code{system.file("inst",package="isobar")} contains R, Sweave, and LaTeX
  files as examples of how to create XLS and PDF reports using \code{isobar}.
  \describe{
    \item{create_reports.R}{Call with Rscript. It is the main file which 
    1) calls a perl script to generate
    a XLS report, 2) generates a LaTeX quality control and analysis report.
    for the XLS report the script pl/tab2xls.pl is used, which concetenates 
    CSV files to a XLS. See Perl requirements. Sweave is called on 
    report/isobar-qc.Rnw and report/isobar-analysis.Rnw. All files are written
    the working directory.}
    \item{isobar-qc.Rnw}{Quality control Sweave file.}
    \item{isobar-analysis.Rnw}{Data analysis Sweave file.}
    \item{properties.conf}{Default configuration for data analysis.}
    \item{report-utils.R}{Helper functions which are used in the Sweave files
    and generate TeX output}
    \item{report-utils.tex}{LaTeX functions for plotting tikz graphics, etc.}
  }
}

\author{Florian P. Breitwieser, Jacques Colinge}

\seealso{
  \link{IBSpectra},
  \link{isobar-preprocessing}
  \link{isobar-analysis}
}

