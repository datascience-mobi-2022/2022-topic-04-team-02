\name{mouse4302mmenstprobe}
\docType{data}
\alias{mouse4302mmenstprobe}
\title{Probe sequence for microarrays of type mouse4302mmenst.}
\description{This data object was automatically created by the package AnnotationForge version 1.32.0.}
\usage{data(mouse4302mmenstprobe)}
\format{A data frame with 754530 rows and 6 columns, as follows.
  \tabular{lll}{
    sequence \tab character \tab probe sequence \cr
    x  \tab integer  \tab x-coordinate on the array \cr
    y  \tab integer  \tab y-coordinate on the array \cr
    Probe.Set.Name      \tab character \tab Affymetrix Probe Set Name \cr
    Probe.Interrogation.Position \tab integer \tab Probe Interrogation Position \cr
    Target.Strandedness \tab factor \tab Target Strandedness\cr
}}

\source{The probe sequence data was obtained from http://www.affymetrix.com.}

\keyword{datasets}
\examples{
  mouse4302mmenstprobe
  as.data.frame(mouse4302mmenstprobe[1:3,])
}
