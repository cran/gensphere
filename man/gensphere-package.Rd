\name{gensphere-package}
\alias{gensphere-package}
\docType{package}
\title{gensphere
}
\description{Define and compute with generalized spherical distributions - multivariate probability
 laws that are specified by a star shaped contour (directional behavior) and a radial component.
}
\details{
This package implements some classes of generalized spherical distributions in dimensions 2, 3, and above.
Functions \code{cfunc.new}, \code{cfunc.add.term}, \code{cfunc.finish} give a flexible way to define a range of shapes for the 
star-shaped contours.  Then function \code{gensphere} defines a generalized spherical distribution
using a contour function and a specification of the radial term.  Function \code{dgensphere} is used 
to compute the multivariate density g(x) for X and function  \code{rgensphere} is
used to simulate a sample random vectors with the (approximate) distribution X.

A large class of distribution can be described as generalized spherical laws. 
In particular, all isotropic/radially symmetric distributions and all elliptically contoured
distributions are generalized spherical laws. Such distributions can be represented as:
\eqn{\mathbf{X} = R \mathbf{S}}{X = R S,}
where R is a positive random variable and \bold{S} is a random vector distributed uniformly (with respect to surface area) on 
the contour, see Nolan (2015).  

Throughout this package, points in d-dimensional space are represented as column vectors; this is different
than what base \R and packages \pkg{mvmesh}, \pkg{geometry}, etc. use; but it is the same as package \pkg{SphericalCubature},
\pkg{SimplicialCubature}, and other packages.

This research was supported by an agreement with Cornell University, Operations 
Research & Information Engineering, under contract W911NF-12-1-0385 from the Army 
Research Development and Engineering Command.


Please let me know if you find any mistakes.  I will try to fix bugs promptly.
Constructive comments for improvements are welcome;  
actually implementing any suggestions will be dependent on time constraints.

}

\references{
B. C. Arnold, E. Castillo and J. M. Sarabia, Multivariate distributions defined in terms of contours,
J. Stat. Planning and Inference, 138, 4158 - 4171, 2008

C. Fernandez, J. Osiewalksi and M. F. J. Steel, Modeling and Inference with v-Spherical Distributions,
J. Amer. Stat. Assoc., 90, 1331-1340, 1995
    
J. P. Nolan, Models for generalized spherical and related distributions.
\href{http://arxiv.org/abs/1510.06772}{arXiv preprint}, Sept. 2015

}
\author{
John P Nolan

Maintainer: John P Nolan

}
\seealso{
\code{\link{cfunc.new}}, \code{\link{gensphere}} 
}