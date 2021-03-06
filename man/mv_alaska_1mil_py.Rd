% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TrendNPS.r
\docType{data}
\name{mv_alaska_1mil_py}
\alias{mv_alaska_1mil_py}
\title{Alaska Coastline 1:1,000,000 - Digital Chart of the World (DCW).}
\format{A shapefile with 811 constituent polygons. To use, run

\code{AK <- rgdal::readOGR(paste0(find.package("TrendNPS"),"/extdata"),"mv_alaska_1mil_py")}.

\describe{The shapefile contains the polygonal extent of the US state of Alaska.  This coverage was taken 
from ESRI's Digital Chart of the World and clipped to the State of Alaska. The data was extracted from 
the PONET coverage from each tile. For more information, see the DCW Data Dictionary.   

The attribute table contains the following variables, whose rows 
contain the same data:

\itemize{ 
  \item TYPE. Uniformly "A" throughout.  Unknown meaning. 
  \item REFRESH_01. Uniformly "2011-03-23" throughout.  Unknown meaning.  
  \item POLNSTAT. Uniformly 0.0000 throughout.  Unknown meaning.  
  \item POLNTYPE. Uniformly 0.0000 throughout.  Unknown meaning.  
}
}}
\description{
Alaska Coastline 1:1,000,000 - Digital Chart of the World (DCW).
}
\references{
Alaska Department of Natural Resources, Information Resource Management.

Accessed January 13, 2017 and available at \url{http://asgdc.alaska.gov#2972}.
}
\author{
Leigh Ann Starcevich of Western EcoSystems Technology, Inc.
}
