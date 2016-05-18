\name{cfunc.new}
\alias{cfunc.new}
\alias{cfunc.add.term}
\alias{cfunc.finish}
\alias{cfunc.eval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Define and evaluate a contour function
}
\description{
The directional part of a generalized spherical distribution is defined
by a contour function, cfunc for short.  These functions are used to
define a contour function and then evaluate it.
}
\usage{
cfunc.new(d)
cfunc.add.term(cfunc, type, k)
cfunc.finish(cfunc, nsubdiv = 2, maxEvals=100000, ...)
cfunc.eval(cfunc, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{dimension of the space
}
  \item{cfunc}{an object of class "gensphere.contour"
}  
  \item{type}{string describing what type of term to add to the contour
}
  \item{k}{a vector of constants used to specify a term
}  
  \item{x}{a (d x n) matrix, with columns x[,i] being points in R^d
}    
  \item{nsubdiv}{number of dyadic subdivisions of the sphere, controls the refinement of the tessellation. Typical values are 2 or 3.
}
  \item{maxEvals}{maximum number of evaluations of the integrand function, see details section below
}
  \item{...}{optional arguments to pass to integration routine, see details section below
}

}
\details{A contour function describes the directional behavior of a generalized spherical distribution.
In this package, a contour function is built by calling three functions:
\code{cfunc.new} to start the definition of a d-dimensional contour, \code{cfunc.add.term} to add a
new term to a contour (may be called more than once), and \code{cfunc.finish} to complete the
defintion of a contour.  

When adding a term, type is one of the literal strings "constant", "elliptical", "proj.normal", "lp.norm", "gen.lp.norm", or "cone".
The vector k contains the constants necessary to specify the desired shape.  k[1]=the first element of k
is always a scale, it allows one to expand or contract a 
shape.  The remaining elements of k depend on the value of 'type'.:

\itemize{
\item{"constant": }{k is a single number; k[1]=radius of the sphere}
\item{"elliptical": }{k is a vector of length d^2+1; k[1]=scale, k[2:(d^2+1)] specify the symmetric positive definite matrix B with is used to compute
the elliptical contour}
\item{"proj.normal": }{k is a vector of length d+2; k[1]=scale, mu=k[2:(d+1)] is the vector pointing in the direction where the normal bump is centered,
    k[d+2]=scale of the isotropic normal bump}
\item{"lp.norm": }{k is a vector of length 2; k[1]=scale and k[2]=p, the power used in the l_p norm}
\item{"gen.lp.norm": }{k is vector of length 2+m*d for some positive integer m; k[1]=scale, k[2]=p, the power used in the l-p norm, k[3:(2+m*d)] 
  spcifies a matrix A that is used to compute || A x ||_p }
\item{"cone": }{k is a vector of length d+2, k[1]=scale, mu=k[2:(d+1)]= the center of the cone, k[d+2]=base of the cone}
}

Note that \code{cfunc.finish} does a lot of calculation, and may take a while, especially
in dimension d > 2.  The most time consuming part is numerically integrating over the contour,
a (d-1) dimensional surface in d-dimensional space and in tesselating the contour in a way that
focuses on the bulges in the contour from cones and normal bumps.  The integration is 
required to calculate the norming constant needed to compute the density.  This integration is
performed by using function \code{adaptIntegrateSphereTri} in \pkg{SphericalCubature} and is 
numerically challenging.  In dimension d > 3 or if nsubdiv > 4, users may need to adjust the arguments maxEvals
and ...  The default value maxEvals=100000 workw in most 3 dim. problems, and it takes a
few seconds to calculate.  (For an idea of the size and time required, a d=4 dim. case
used maxEvals=1e+7 and took around 5 minutes. A d=5 dim. case used maxEvals=1e+8, used
160167 simplices and took over 2 days.)  Note that this calculation is only done once;
calculating densities and simulating is still fast in higher dimensions.  It may be useful
to save a complicated/large contour object so that it can be reused across
R sessions via \code{save(cfunc)} and \code{load(cfunc)}.

Note: the first time \code{cfunc.finish} is called, a warning message 
about "no degenerate regions are returned" is printed by the package
\pkg{geometry}.  I do not know how to turn that off; so just ignore it.

\code{cfunc.eval} is used to evaluate a completed
contour function.
}
\value{\code{cfunc.new} and \code{cfunc.add.term} return a list that is an incomplete definition
of a contour function.  \code{cfunc.finish} completes the definition and returns an S3 object of
class "gensphere.contour" with fields:

\item{d}{dimension}
\item{m}{number of terms in the contour function, i.e. the number of times \code{cfunc.add.term} was called}
\item{term}{a vector length m of type list, with each list describing a term in the contour function}
\item{norm.const}{norming constant}
\item{functionEvaluations}{number of times the integrand (contour) function is evaluated by \code{adaptIntegrateSphereTri}}
\item{tessellation}{an object of type "mvmesh" that gives a geometrical description of the contour. It is used
 to plot the contour and to simulate from the contour}
\item{tessellation.weights}{weights used in simulation; surface area of the simplices in the tessellation}
\item{simplex.count}{vector of length 3 giving the number of simplices used at the end of
  three internal stages of construction: after initial subdivision, after refining the sphere based on cones and bumps, 
    and final count have adaptive integration routine partitions the sphere}
\code{cfunc.eval} returns a vector of length n=nrow(x); y[i] = cfunc(x[,i]) = value of the contour function at point x[,i].

The plots below show the three contours defined in the examples below.

\if{html}{
\figure{contours2d.png}

\figure{star.png}
}
\if{latex}{
\figure{contours2d.png}{options: width=4in}

\figure{star.png}{options: width=3in}
}

}
\examples{
# 2-dim diamond
cfunc1 <- cfunc.new(d=2)
cfunc1 <- cfunc.add.term( cfunc1,"lp.norm",k=c(1,1))
cfunc1 <- cfunc.finish( cfunc1 )
cfunc1
cfunc.eval( cfunc1, c(sqrt(2)/2, sqrt(2)/2) )
\donttest{plot( cfunc1, col='red', lwd=3, main="diamond contour")}

# 2-dim blob
cfunc2 <- cfunc.new(d=2)
cfunc2 <- cfunc.add.term( cfunc2,"constant",k=1)
cfunc2 <- cfunc.add.term( cfunc2,"proj.normal",k=c( 1, sqrt(2)/2, sqrt(2)/2, 0.1) )
cfunc2 <- cfunc.add.term( cfunc2,"proj.normal",k=c( 1, -1,0, 0.1) )
cfunc2 <- cfunc.finish( cfunc2, nsubdiv=4 )
\donttest{
plot( cfunc2, col='green', lwd=3, main="contour with two bumps")

# 3-dim star-like contour
cfunc3 <- cfunc.new(d=3)
cfunc3 <- cfunc.add.term( cfunc3, "constant", k=0.5 )
for (i in 1:3) {
  u <- c(0,0,0);  u[i] <- 1  # i-th unit vector
  cfunc3 <- cfunc.add.term( cfunc3, "cone", k=c(1,u,.2) )
  cfunc3 <- cfunc.add.term( cfunc3, "cone", k=c(1,-u,.2) )
}  
cfunc3 <- cfunc.finish( cfunc3, maxEvals=400000 ) # this is slow, ~ 15 min

plot( cfunc3, show.faces=TRUE, col='blue')
nS <- dim(cfunc3$tessellation$S)[3]
title3d( paste("star with",nS,"simplices"))
}
}
