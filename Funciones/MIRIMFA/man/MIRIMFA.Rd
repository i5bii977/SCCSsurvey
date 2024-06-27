\name{MIRIMFA}
\alias{MIRIMFA}

\title{
Multiple Imputation for Regularised Iterative Multiple Factor Analysis (MIRIMFA)
}

\description{
This function performs Multiple Imputation by Regularised Iterative Multiple Factor Analysis (MIRIMFA) on a dataset with missing values. It imputes missing values using the Regularised Iterative Multiple Factor Analysis (MFA) (Husson & Josse, 2013) method and creates multiple imputations to track the uncertainty of the imputation. The imputed dataset is then analysed using the MFA method to compute a compromise space in the fashion of Viollet et al. (2016).
}

\usage{
MIRIMFA(X,
        group,
        ncp = 2,
        type = rep("s", length(group)),
        M = 5,
        row.w = NULL,
        coeff.ridge = 1,
        threshold = 1e-06,
        ind.sup = NULL,
        num.group.sup = NULL,
        maxiter = 1000,
        name.group = NULL)
}

\arguments{
  \item{X}{Data matrix with missing values.}
  \item{group}{Vector specifying the group membership of the variables in the analysis.}
  \item{ncp}{Number of dimensions for Multiple Factor Analysis.}
  \item{type}{Type of analysis for each group of variables.}
  \item{M}{Number of multiple imputations to generate.}
  \item{row.w}{Weights for rows.}
  \item{coeff.ridge}{Ridge coefficient for Regularised MFA.}
  \item{threshold}{Threshold for stopping criterion.}
  \item{ind.sup}{Indices of supplementary individuals.}
  \item{num.group.sup}{Number of supplementary groups.}
  \item{maxiter}{Maximum number of iterations for imputation.}
  \item{name.group}{Optional. Names for each group of variables.}
}

\details{
The function iteratively performs the following steps for each multiple imputation:

1. Imputes missing values using the Regularized MFA method.
2. Computes the compromise space using the STATIS method.
3. Constructs a total disjunctive table for the imputed dataset.
4. Imputes missing values in the total disjunctive table.
5. Returns the multiple imputations, compromise space, and other relevant results.
}

\references{
  F. Husson, J. Josse (2013) Handling missing values in multiple factor analysis. Food Quality and Preferences, 30 (2), 77-85.
Josse, J. and Husson, F. missMDA (2016). A Package for Handling Missing Values in Multivariate Data Analysis. Journal of Statistical Software, 70 (1), pp 1-31 tools:::Rd_expr_doi("10.18637/jss.v070.i01").
  Voillet V., Besse P., Liaubet L., San Cristobal M., González I. (2016). Handling missing rows in
multi-omics data integration: Multiple Imputation in Multiple Factor Analysis framework. BMC
Bioinformatics, 17(40).
}

\author{
Rodrigo Olivera (rmmolivera@tutanota.com)


\seealso{
\code{\link{plotFactorsMIRIMFA}},
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
