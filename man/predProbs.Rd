\name{predProbs}
\alias{predProbs}
\title{Predicted Probabilities for a Variety of Logit and Probit Models}
\description{This an R function for computing predicted probabilities
for binary & ordinal logit and probit, (partial) generalized ordinal & multinomial   logit models estimated with glm, clm (in \pkg{ordinal}), and vglm (in \pkg{VGAM}) commands. It returns a data frame with each column containing the predicted probabilities for a specific response y value given a set of chosen independent variable
settings.}
\usage{
predProbs(model, specs, method="logit")
}
\arguments{
\item{model}{An input model object estimated with glm, clm, or vglm. The order
of the right-hand variables must be the same as that of the "specs" argument.
}

\item{specs}{A data frame with each row specifying the chosen values of
all the independent variables in the model.
}

\item{method}{Default "logit"; altenative methods are "probit," "mlogit," and "gologit".The "logit" and "probit" method can be estimated with \pkg{glm} or \pkg{clm} of the \pkg{ordinal} package while "mlogit" and "gologit" can be estimated with \pkg{vglm} in the \pkg{VGAM} package. For multinomial logit models, use the last choice as the reference category.
}
}

\value{
The function outputs a data frame of J number of columns, with each column containing the predicted probabilities p(y=j) with j = 1, ..., J for ordinal models, j = 1, 0 for binary models, and j = 1, ..., Ref for multinomial models. The rows are defined the same as in the input "specs" argument.
}

\examples{
df <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")
df$hRank[df$rank==1 | df$rank==2] <- 1
df$hRank[df$rank==3 | df$rank==4] <- 0
logit1 <- glm(admit ~ gre + hRank + gpa, data=df, binomial)
setval1 <- expand.grid(gre=seq(250,800,50), hRank=0:1, gpa=mean(df$gpa))
predprobs1 <- predProbs(logit1, setval1)
probit1 <- glm(admit ~ gre + hRank + gpa, data=df, binomial(link=probit))
predprobs2 <- predProbs(probit1, setval1, method="probit")

}

\references{
Tim F. Liao, 1994. \emph{Interpreting Probability Models: Logit, Probit, and Other Generalized Linear Models}. Thousand Oaks, CA: Sage.

J. Scott Long, 1997. \emph{Regression Models for Categorical and Limited Dependent Variables}. Thousand Oaks, CA: Sage.
}
