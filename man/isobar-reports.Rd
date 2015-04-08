\name{isobar-reports}

\alias{isobar-reports}
\alias{connect.nodes}
\alias{create.reports}
\alias{create.meta.reports}
\alias{draw.boxplot}
\alias{draw.protein.group}
\alias{draw.proteingroup.row}
\alias{initialize.env}
\alias{load.properties}
\alias{modifs}
\alias{print_longtablehdr}
\alias{tikz.proteingroup}
\alias{write.xls.report}
\alias{print_longtablehdr_peptide}
\alias{transform_pepmodif}
\alias{property}
\alias{write.tex.commands}
\alias{load.tex.properties}

\alias{get_n_proteins}
\alias{print_classlabels_tbl}
\alias{print_groupsize}
\alias{print_protein_grp_info}
\alias{print_protein_grp_tbl}
\alias{print_protein_notquant_tbl}
\alias{print_protein_quant_tbl}
\alias{print_sign_proteins_tbl}

\alias{testPdflatex}
\alias{testPerl}

\title{Isobar reports}

\description{
  Generation of LaTeX and XLS reports is helped with
  functions which facilitate the gathering of relevant information and
  creation of tikz plots. \code{create.reports} parses properties (by
  calling \code{load.properties}) and initialize environments and
  computations (by calling \code{initialize.env}) required by the
  reports, calls Sweave and pdflatex.
}

\usage{
create.reports(properties.file = "properties.R", 
               global.properties.file = system.file("report","properties.R", package = "isobar"),
               args = NULL, ...,
               recreate.properties.env = TRUE, recreate.report.env = TRUE)

load.properties(properties.file = "properties.R",
                global.properties.file = system.file("report","properties.R",package="isobar"),
                args = NULL, ...)

initialize.env(env, properties.env)

}

\arguments{
  \item{properties.file}{File which holds the parameters for data
    analysis and report generation. It is parsed as R code after the
    global report configuration file \code{global.properties.file}  
    and defines peaklists, identification files, significance levels,
    etc. See the global properties file for the available options and
    values.}
  \item{global.properties.file}{
    \code{system.file("report","properties.R",package="isobar")}}
  \item{args}{Additional (command line) arguments which overrids those
    in \code{properties.file}. }
  \item{...}{Additional properties.}
  \item{recreate.properties.env}{Whether a properties.env existing in the global
  environment should be used, or it should be recreated.}
  \item{recreate.report.env}{Whether a report.env existing in the global
  environment should be used, or it should be recreated.}
  \item{env}{Item to be initialized.}
  \item{properties.env}{Environment into which properties are read.}
}

\details{
  The directory \code{inst} in the isobar installation directory 
  \code{system.file("inst",package="isobar")} contains R, Sweave, and LaTeX
  files as examples of how to create XLS and PDF reports using \code{isobar}.
  \describe{
    \item{create_reports.R}{Call with Rscript. It is the main file which
      \enumerate{
	\item{parses command line options. \code{--compile} and
	  \code{--zip} are parsed directly and given as arguments to
	  \code{create.reports}. Other arguments are given
	  \code{\link{load.properties}}.}	
	\item{calls a perl script to generate a XLS report}
	\item{generates a LaTeX quality control and analysis report}
	}
	for the XLS report the script pl/tab2xls.pl is used, which
	concetenates CSV files to a XLS. See Perl requirements.
	Sweave is called on report/isobar-qc.Rnw and
	report/isobar-analysis.Rnw. All files are written the
	working directory.}
    \item{isobar-qc.Rnw}{Quality control Sweave file.}
    \item{isobar-analysis.Rnw}{Data analysis Sweave file.}
    \item{properties.R}{Default configuration for data analysis.}
    \item{report-utils.tex}{LaTeX functions for plotting tikz graphics, etc.}
  }
}

\author{Florian P Breitwieser}

\seealso{
  \link{IBSpectra},
  \link{isobar-preprocessing}
  \link{isobar-analysis}
}

