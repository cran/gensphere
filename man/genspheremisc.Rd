\name{genspheremisc}
\alias{gs.cone}
\alias{gs.elliptical}
\alias{gs.gen.lp.norm}
\alias{gs.lp.norm}
\alias{gs.pdf2d.plot}
\alias{gs.proj.normal}
\alias{gs.vfunc.eval}
\alias{RefineSphericalTessellation}
\alias{NearbyPointsOnSphere}
\alias{RotateInNDimensions}
\alias{print.gensphere.contour}
\alias{plot.gensphere.contour}
\alias{print.gensphere.distribution}

\title{Miscellaneous functions used in working with generalized spherical laws
}
\description{Miscellaneous internal functions for the gensphere package.
}
\usage{
gs.cone(x, mu, theta0)
gs.elliptical(x, B)
gs.gen.lp.norm(x, p, A)
gs.lp.norm(x,p)
gs.proj.normal(x, mu, sigma)
gs.vfunc.eval(cfunc, x)
gs.pdf2d.plot(gs.dist, xy.grid = seq(-10, 10, 0.1) )
RefineSphericalTessellation(V1, V2)
NearbyPointsOnSphere(x, epsilon)
RotateInNDimensions(x, y)
\method{print}{gensphere.contour}(x,...)
\method{print}{gensphere.distribution}(x,...)
\method{plot}{gensphere.contour}(x,multiplier=1,...)
}
\arguments{
  \item{x,y}{vectors representing points in d-dimensional space
}
  \item{mu}{direction of the mode for a cone/normal bump
}
  \item{theta0}{angle between peak of the cone and the base of the cone
 }
  \item{B}{(d x d) positive definite shape matrix
}
  \item{A}{ matrix used to compute || A x ||_p
}
  \item{p}{power of the l^p norm; p=2 is Euclidean distance
}
  \item{gs.dist}{object of class "gensphere.distribution" defined by \code{gensphere}
}
  \item{xy.grid}{a matrix of (x,y) values in 2-dimensions
}
  \item{cfunc}{an object of class "gensphere.contour" defined by \code{cfunc.new}, etc.
}
  \item{\dots}{optional arguments to the 2-dimensional plot, e.g. col='red', etc.
}
  \item{sigma}{scale parameter for a normal bump
}
  \item{epsilon}{vector of positive numbers where there are points added around a particular direction
}
  \item{V1,V2}{matrices of vertices which are joined together to get a refinement of the grid
}
  \item{multiplier}{a positive number used to scale the contour
}
}

\details{These are undocumented functions that are used internally.  The functions 
\code{gs.cone}, \code{gs.elliptical}, \code{gs.gen.lp.norm}, \code{gs.lp.norm}, 
\code{gs.proj.normal}, \code{gs.vfunc.eval} are used in evalutating a contour function.
\code{RefineSphericalTessellation}, \code{NearbyPointsOnSphere} are used in defining 
the tessellation of the contour that identifies bumps and cones.
\code{gs.pdf2d.plot} and the plot/print methods are initial attempts at plotting and
printing a summary of objects. 

These functions may change or disappear in the future.
}

