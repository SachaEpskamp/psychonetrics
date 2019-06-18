\name{StarWars}
\alias{StarWars}
\docType{data}
\title{
Star Wars dataset
}
\description{
This questionaire was constructed by Carolin Katzera, Charlotte Tanis, Esther Niehoff, Myrthe Veenman, and Jason Nak as part of an assignment for a course on confirmatory factor analysis (http://sachaepskamp.com/SEM2018). They also collected the data among fellow psychology students as well as through social media. 
}
\usage{data("StarWars")}
\format{
  A data frame with 271 observations on the following 13 variables.
  \describe{
    \item{\code{Q1}}{I am a huge Star Wars fan! (star what?)}
    \item{\code{Q2}}{I would trust this person with my democracy <picture of Jar Jar Binks>}
    \item{\code{Q3}}{I enjoyed the story of Anakin's early life}
    \item{\code{Q4}}{The special effects in this scene are awful <video of the Battle of Geonosis>}
    \item{\code{Q5}}{I would trust this person with my life <picture of Han Solo>}
    \item{\code{Q6}}{I found Darth Vader'ss big reveal in "Empire" one of the greatest moments in movie history}
    \item{\code{Q7}}{The special effects in this scene are amazing <video of the Death Star explosion>}
    \item{\code{Q8}}{If possible, I would definitely buy this droid <picture of BB-8>}
    \item{\code{Q9}}{The story in the Star Wars sequels is an improvement to the previous movies}
    \item{\code{Q10}}{The special effects in this scene are marvellous <video of the Starkiller Base firing>}
    \item{\code{Q11}}{What is your gender?}
    \item{\code{Q12}}{How old are you?}
    \item{\code{Q13}}{Have you seen any of the Star Wars movies?}
  }
}
\details{
The questionaire is online at https://github.com/SachaEpskamp/SEM-code-examples/blob/master/CFA_fit_examples/StarWars_questionaire.pdf. The authors of the questionaire defined a measurement model before collecting data: Q2 - Q4 are expected to load on a "prequel" factor, Q5 - Q7 are expected to load in a "originals" factor, and Q8 - Q10 are expected to load on a "sequal" factor. Finally, Q1 is expected to load on all three.
}
\source{
https://github.com/SachaEpskamp/SEM-code-examples/blob/master/CFA_fit_examples
}
\examples{
data(StarWars)
}
\keyword{datasets}
