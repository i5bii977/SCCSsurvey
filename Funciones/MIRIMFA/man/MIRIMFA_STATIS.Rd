\name{MIRIMFA_STATIS}
\alias{MIRIMFA_STATIS}

\title{
Pooling of manual Multiple Imputation by Regularised Iterative Multiple Factor Analysis
}

\description{
The MIRIMFA_STATIS function allows for pooling of multiple imputations using the Regularised Iterative Multiple Factor Analysis (RIMFA) method through the function MIRIMFAextraimp. It is particularly useful for large databases with many missing values where running the MIRIMFA function directly might be resource-intensive. This function enables users to take the imputations saved from MIRIMFAextraimp and pool the imputed datasets to obtain the same results as MIRIMFA.
}

\usage{
MIRIMFAextraimp(X,
folder_path,
prefix_imputations = "RIMFA_number_")
}

\arguments{
\item{X}{A data frame containing the variables to be imputed.}
\item{folder_path}{A character string specifying the path to the imputed datasets.}
\item{prefix_imputations}{Optional. A character string specifying the prefix of the objects to be taken from the path.}
}

\details{
The MIRIMFA_STATIS function pools multiple imputations using the Regularised Iterative Multiple Factor Analysis (RIMFA) method with iterative imputation control. It takes the imputed dataset from a folder or directory, and calculates a compromise space to obtain the final results of the multiple imputation in a similar way to MIRIMFA. The benefit of using MIMIMFAextraimp together with MIRIMFA_STATIS is that one can interrupt the process and retake at any moment, which is useful for large datasets with a big proportion of missing values to be able to rationalise computational resources.
}

\author{
Rodrigo Olivera (rmmolivera@tutanota.com)
}

\seealso{
\code{\link{MIRIMFA}},
\code{\link{MIRIMFAextraimp}}
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
