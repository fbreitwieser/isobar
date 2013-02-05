\name{isobar.data}
\alias{ibspiked_set1}
\alias{ibspiked_set2}
\alias{noise.model.hcd}
\alias{ib_phospho}
\docType{data}
\title{
  Isobar Data packages
}
\description{
  ibspiked_set1 and ibspiked_set2 are objects of class iTRAQ4plexSpectra.
  It contains over 160 protein groups, over 1600 peptides from about 15,000 spectra each,
  mainly from background proteins and three spiked-in Ceruplasmins 
  (CERU_HUMAN, CERU_MOUSE, CERU_RAT).
}
\usage{
  data(ibspiked_set1)
  data(ibspiked_set2)
  data(ib_phospho)
}
\format{
  iTRAQ4plexSpectra objects.
}
\source{
  isobar publication. Acquired on Orbitrap instrument w/ 20 offline-fractions 
  and HCD fragmentation.
}
\examples{
data(ibspiked_set1)
print(ibspiked_set1)
}
\keyword{datasets}
