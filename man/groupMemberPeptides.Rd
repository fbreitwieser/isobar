\name{groupMemberPeptides}
\alias{groupMemberPeptides}
\title{
    Peptide info for protein group members
}
\description{
    For a given reporter protein group identifier, information on
    its peptides is returned. It contains information on how the
    peptides are shared and in which member they occur.
}
\usage{
groupMemberPeptides(x, reporter.protein.g, ordered.by.pos = TRUE, only.first.pos = TRUE)
}
\arguments{
  \item{x}{ProteinGroup object}
  \item{reporter.protein.g}{group reporter protein}
  \item{ordered.by.pos}{if TRUE, start position of peptides in 
      proteins is exported and peptides are ordered by position
}
  \item{only.first.pos}{
      if TRUE, only first occurence of peptide in protein is reported
}
}
\value{
    list of two:
        [1] peptide.info: data.frame
             peptide
             specificity
             n.shared.groups
             n.shared.proteins
             start.pos
        [2] group.member.peptides: data.frame
             each column corresponds to a group member, and
             each row to a peptide
}
\author{Florian P Breitwieser}

\examples{
    data(ibspiked_set1)
    protein.group <- proteinGroup(ibspiked_set1)
    ceru.rat <- protein.g(protein.group,"CERU_RAT")
    groupMemberPeptides(protein.group,ceru.rat)

## find protein groups with members
    t <- table(proteinGroupTable(protein.group)$reporter.protein)
    t[t>2]
    protein.g <- names(t)[t>2][1]
    groupMemberPeptides(protein.group,protein.g)
}
