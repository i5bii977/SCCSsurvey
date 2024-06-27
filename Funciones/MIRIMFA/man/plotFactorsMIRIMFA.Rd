\name{plotFactorsMIRIMFA}
\alias{plotFactorsMIRIMFA}

\title{
Plot Multiple Imputation Results by Regularised Iterative Multiple Factor Analysis (MIRIMFA)
}

\description{
This function plots the results of Multiple Imputation by Regularised Iterative Multiple Factor Analysis (MIRIMFA). It visualizes the projection of the main factors of each imputed dataset onto the dimensions of the reference compromise space.
}

\usage{
plotFactorsMIRIMFA(object,
                    comp = c(1, 2))
}

\arguments{
  \item{object}{An object of class 'MIRIMFA', obtained from the MIRIMFA function.}
  \item{comp}{A vector specifying the components to plot (default is c(1,2)).}
}

\details{
The plotFactorsMIRIMFA function generates plots to visualize the results of Multiple Imputation by Regularised Iterative Multiple Factor Analysis (MIRIMFA). It displays the projection of the main factors of each imputed dataset onto the dimensions of the reference compromise space.

The main purpose of this plot is to represent the multiple imputed datasets on a reference compromise configuration, which is obtained from the MIRIMFA analysis conducted on the incomplete dataset.
}

\references{
  Audigier, V., Husson, F., Josse, J. (2016). MIMCA: Multiple imputation for categorical variables with multiple correspondence analysis
  F. Husson, J. Josse (2013) Handling missing values in multiple factor analysis. Food Quality and Preferences, 30 (2), 77-85.
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
\code{\link{plotVarMIRIMFA}},
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
