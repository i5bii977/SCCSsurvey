\name{plotVarMIRIMFA}
\alias{plotVarMIRIMFA}

\title{
Plot Variables in Multiple Imputation Results by Regularised Iterative Multiple Factor Analysis (MIRIMFA)
}

\description{
This function generates scatter plots to visualize the variables in the results of Multiple Imputation by Regularised Iterative Multiple Factor Analysis (MIRIMFA). It displays the relationships between variables and the projection of each variable onto the dimensions of the reference compromise space.
}

\usage{
plotVarMIRIMFA(object,
comp = c(1, 2))
}

\arguments{
\item{object}{An object of class 'MIRIMFA', obtained from the MIRIMFA function.}
\item{comp}{A vector specifying the components to plot (default is c(1,2)).}
}

\details{
The plotVarMIRIMFA function generates scatter plots to visualize the relationships between variables in the results of Multiple Imputation by Regularised Iterative Multiple Factor Analysis (MIRIMFA). It displays the projection of each variable from an MFA of each imputation onto the dimensions of the reference compromise space.
}

\references{
  Josse, J. and Husson, F. missMDA (2016). A Package for Handling Missing Values in Multivariate Data Analysis. Journal of Statistical Software, 70 (1), pp 1-31 tools:::Rd_expr_doi("10.18637/jss.v070.i01").
  Voillet V., Besse P., Liaubet L., San Cristobal M., González I. (2016). Handling missing rows in
multi-omics data integration: Multiple Imputation in Multiple Factor Analysis framework. BMC
Bioinformatics, 17(40).
}

\author{
Rodrigo Olivera (rmmolivera@tutanota.com)
}

\seealso{
\code{\link{MIRIMFA}},
\code{\link{plotFactorsMIRIMFA}},
\code{\link{MIRIMFAextraimp}},
\code{\link{MIRIMFA_STATIS}}
}


\examples{
# Example usage:
data (poison)
# Amputation of data
total_cells <- length(as.matrix(poison))
random_rows <- sample(1:nrow(poison), round(0.2 * total_cells), replace = TRUE)
random_cols <- sample(1:ncol(poison), round(0.2 * total_cells), replace = TRUE)
poison_amp <- poison
for (i in 1:length(random_rows)) {
  poison_amp[random_rows[i], random_cols[i]] <- NA
}
# Running MIRIMFA
MIpoison <- MIRIMFA(poison_amp, group=c(2,2,5,6), type=c("s","n","n","n"), M = 20)

# Plotting variables and Factors
plotVarMIRIMFA(MIpoison)
plotFactorsMIRIMFA(MIpoison)
}
