\name{MIRIMFAextraimp}
\alias{MIRIMFAextraimp}

\title{
Multiple Imputation by Regularised Iterative Multiple Factor Analysis with Iterative Imputation extra imputations
}

\description{
The MIRIMFAextraimp function allows for multiple imputations using the Regularised Iterative Multiple Factor Analysis (RIMFA) method with iterative imputation control. It is particularly useful for large databases with many missing values where running the MIRIMFA function directly might be resource-intensive. This function enables users to perform imputation one iteration at a time, save each imputation, and then pool the imputed datasets later using MIRIMFA_STATIS to obtain the same results as MIRIMFA.
}

\usage{
MIRIMFAextraimp(X,
group,
ncp = 2,
type = rep("s", length(group)),
manual_threshold = 1,
row.w = NULL,
coeff.ridge = 1,
threshold = 1e-06,
ind.sup = NULL,
num.group.sup = NULL,
maxiter = 1000,
folder_path = here("imputationsRIMFA"),
name.group = NULL,
start = 0,
initial_seed = 1)
}

\arguments{
\item{X}{A data frame containing the variables to be imputed.}
\item{group}{A vector specifying the grouping structure of the variables.}
\item{ncp}{An integer specifying the number of components to use in the analysis.}
\item{type}{A character vector specifying the types of variables in each group.}
\item{manual_threshold}{An integer specifying the number of imputations to perform before stopping.}
\item{row.w}{Optional. A matrix of weights for the rows of \code{X}.}
\item{coeff.ridge}{A numeric value specifying the ridge coefficient to use in the imputation.}
\item{threshold}{A numeric value specifying the convergence threshold for the iterative imputation process.}
\item{ind.sup}{Optional. A list specifying supplementary individuals.}
\item{num.group.sup}{Optional. A list specifying the number of supplementary groups.}
\item{maxiter}{An integer specifying the maximum number of iterations for the imputation process.}
\item{folder_path}{A character string specifying the path to save the imputed datasets.}
\item{name.group}{Optional. A character vector specifying the names of the groups.}
\item{start}{An integer specifying the starting index for imputations.}
\item{initial_seed}{An integer specifying the initial seed for the random number generator.}
}

\details{
The MIRIMFAextraimp function performs multiple imputations using the Regularised Iterative Multiple Factor Analysis (RIMFA) method with iterative imputation control. It iteratively imputes missing values in the dataset, saves each imputed dataset, and stops once the specified number of imputations is reached. Users can later pool the imputed datasets using MIRIMFA_STATIS to obtain the final results.
}

\author{
Rodrigo Olivera (rmmolivera@tutanota.com)
}

\seealso{
\code{\link{MIRIMFA}},
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
# Running MIRIMFA through manual imputations
start <- 0
initial_seed <- 1
manual_threshold <- 1
for (i in 1:20){
  imputeMIRIMFApoison <- MIRIMFAextraimp(poison_amp,
                          group=c(2,2,5,6),
                          ncp = 7,
                          type=c("s","n","n","n"),
                          manual_threshold = manual_threshold,
                          folder_path = here("RIMFAsPoisonTest"),
                          name.group=c("desc","desc2","symptom","eat"),
                          start = start,
                          initial_seed = initial_seed)
  initial_seed <- imputeMIRIMFApoison$final_seed
  start <- imputeMIRIMFApoison$final_index
  manual_threshold <- manual_threshold + 1
}

MIpoison <- MIRIMFA_STATIS(poison_amp, here("RIMFAsPoisonTest"))
# Plotting variables and Factors
plotVarMIRIMFA(MIpoison)
plotFactorsMIRIMFA(MIpoison)
}
