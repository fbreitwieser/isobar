\name{isobar.data}
\alias{ibspiked_set1}
\alias{noise.model.hcd}
\docType{data}
\title{
  Isobar Data packages
}
\description{
  ibspiked_set1 is a object of class iTRAQ4plexSpectra.
  It contains 161 protein groups, 1653 peptides from nearly 15,000 spectra,
  mainly from background proteins and also three spiked-in Ceruplasmins 
  (CERU_HUMAN, CERU_MOUSE, CERU_RAT).
}
\usage{
  data(ibspiked_set1)
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
