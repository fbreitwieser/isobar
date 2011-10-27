\name{isobar-reports}

\alias{isobar-reports}
\alias{connect.nodes}
\alias{create.reports}
\alias{draw.boxplot}
\alias{draw.protein.group}
\alias{initialize.env}
\alias{load.properties}
\alias{modifs}
\alias{print_longtablehdr}
\alias{tikz.proteingroup}
\alias{write.xls.report}

\title{Isobar reports}

\description{
  Generation of LaTeX and XLS reports is helped with functions
  which facilitate the gathering of relevant information and creation
  of tikz plots. \code{initialize.report} is the main function establishing 
  an environment with all neccessary data and computations.
}

\usage{
create.reports(properties.file = "properties.R", args, report.type = "protein", compile = FALSE)
}

\arguments{
  \item{properties.file}{File which holds the parameters for data analysis and report generation.
    It is written in standard R code and defines peaklists, identification files, significance levels,
    etc. See \code{system.file("report","properties.conf",package="isobar")} for 
    a full set of possible values. It is parsed using \code{\link{sys.source}}.}
  \item{args}{Additional parameters which overrids those in \code{properties.file}. }
  \item{report.type}{Currently, only protein is implemented.}
  \item{compile}{Compile LaTeX source to PDF?}
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
    \item{properties.R}{Default configuration for data analysis.}
    \item{report-utils.tex}{LaTeX functions for plotting tikz graphics, etc.}
  }
}

\author{Florian P. Breitwieser}

\seealso{
  \link{IBSpectra},
  \link{isobar-preprocessing}
  \link{isobar-analysis}
}

