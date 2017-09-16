###################################################################
# R package for generalized spherical distributions
# These are multivariate distributions that generalize spherical
# distributions in that they have level curves that are all scaled
# versions of a common contour.
#         John Nolan, jpnolan@american.edu, March 2013 original version
#          - modified October-November, 2015 to allow simulation from a tessellation 
#            of the contour, turn into a proper R package for CRAN, etc.
#
#####################################################################
cfunc.eval <- function( cfunc, x ) {
#  evaluate r[i] := c( x[,i] ), i=1,...,n, where
#  c(.) is a contour function described in the cfunc object
#  x is a (d x n) matrix, with x[,i] of points in R^d
#
#  This function computes c(x) by combining terms of the
#  following types:
#
# type/name      number of parameters    function in terms of parameters
#    constant          1                c(x) = k[1]
#    elliptical        1+d^2            c(x) = k[1] * sqrt(x^T A x), k[2:(1+d^2)] = A matrix
#    proj.normal       d+2              c(x) = k[1] * g(x,k[2:(d+1)],k[d+2]), where
#                                               g=projected.istropic.normal( )
#    lp.norm           2                c(x) = k[1] * ||x||_p, where p=k[2]
#    gen.lp.norm       2+j*d            c(x) = k[1] * || A x ||_p, where p=k[2] and k[3:(2+j*d)] = A, A a j by d matrix
#    cone              d+2              c(x) = k[1] * (height of a cone with peak 1 at mu=(k[2],k[3],...,k[d+1]))
#                                                and base r=k[d+2]

d <- cfunc$d
if( is.vector(x) ) { x <- matrix( x, nrow=length(x), ncol=1 ) }
stopifnot( nrow(x) == d )
n <- ncol(x)
r <- rep(0.0,n)
rstar <- r
num.star <- 0
for (j in 1:cfunc$m) {
  # evaluate the j-th term in cfunc
  newtype <- charmatch( cfunc$term[[j]]$type, c("constant","elliptical","proj.normal",
           "lp.norm","gen.lp.norm","cone"))
  if (is.na( newtype )) warning("Ignoring unrecognized term in cfunc.eval( )")
  k <- cfunc$term[[j]]$k
  rnew <- switch( newtype,
      rep(1.0,n),
      gs.elliptical( x, matrix(k[-1],d,d,byrow=TRUE) ),
      gs.proj.normal( x,k[2:(d+1)],k[d+2] ),
      gs.lp.norm(x, k[2]),
      gs.gen.lp.norm( x, k[2], matrix( k[-c(1,2)],ncol=d,byrow=TRUE ) ),
      gs.cone( x, matrix( k[2:(d+1)], nrow=d ), k[d+2] )    
  )

  if (cfunc$term[[j]]$invert) {
    rstar <- rstar + k[1]*rnew
    num.star <- num.star + 1
  } else  { r <- r + k[1]*rnew }
}

if (min(rstar) > 0) {
  r <- r + 1/rstar
} else {
  if (num.star > 0) {
    warning("Error in cfunc.eval: some terms have 0 in denominator; 0/0 undefined")
    r <- rep(NA,length(r))
  }
}
return(r) }
#####################################################################

gensphere <- function( cfunc, dradial, rradial, g0 ) {
# define a new generalized spherical distribution
#   cfunc is a contour function defined with cfunc.new, cfunc.add.term, cfunc.finish
#   dradial is a function to simulate from the radial distribution
#   rradial is a function to compute the density of the radial distribution
#   g0 is the value  lim_{r -> 0+} r^(1-d)*dradial(r); if unknown, approximate numerically.
stopifnot( class(cfunc)=="gensphere.contour", is.function(dradial), is.function(rradial), g0 >= 0 )

dist <- list( cfunc=cfunc, dradial=dradial, rradial=rradial, g0=g0 )
class(dist) <- "gensphere.distribution"
return( dist ) }

#####################################################################
cfunc.new <- function( d ) {
# utility function to define a new contour function in dimension d

d <- as.integer( d )
stopifnot( d > 1 )
cfunc <- list( d=d, m=0L, norm.const=NA, term=vector("list",0) )
return( cfunc ) }

#####################################################################
cfunc.add.term <- function( cfunc, type, k) {
# utility function to add a new term to a cfunc list
#  type must be one of the types defined at top
#  k is a vector of doubles used as the parameters to this type of cfunc 

# check that it is a legitimate type and the length of k is correct
d <- cfunc$d
m <- cfunc$m + 1
k <- as.double( k )
invert <- NULL

if (type=="constant") { 
  stopifnot( length(k) == 1 )
  invert <- FALSE
}
if (type=="elliptical") { 
  stopifnot( length(k) == (1+d^2) )
  invert <- TRUE
}
if (type=="proj.normal") { 
  stopifnot( length(k) == (d+2) )
  invert <- FALSE
}
if (type=="lp.norm") { 
  stopifnot( length(k) == 2 )
  invert <- TRUE
}
if (type=="gen.lp.norm") { 
  stopifnot( (length(k)- 2) %% d   == 0 )
  invert <- TRUE
}
if (type=="cone") { 
  stopifnot( length(k) == d+2 )
  invert <- FALSE
}

if (is.null(invert)) stop( "Undefined type in gs.add.cfunc.term" )

cfunc$term[[m]] <- list( type=type,  invert=invert, k=k )
cfunc$m <- m
return(cfunc) }

#####################################################################
cfunc.finish <- function( cfunc, nsubdiv=2, maxEvals=100000, norm.const.method="integrate", ... ) {
# utility function to finish the definition of a contour function

# compute tessellation of contour; start with a sphere
d <- cfunc$d
sphere <- UnitSphere( n=d, k=nsubdiv )

# loop through terms, accumulating new points/vertices to add to the grid
newV <- matrix( 0, nrow=0, ncol=d )
for (i in 1:cfunc$m) {
  cur.type <- (cfunc$term[[i]])$type
  if(  cur.type == "proj.normal" | cur.type == "cone" ) {
    center <- (cfunc$term[[i]])$k[2:(d+1)]
    epsilon <- (cfunc$term[[i]])$k[d+2] * c(0.1,.5,1,1.5,2,2.5)
    a <- NearbyPointsOnSphere( center, epsilon )
    newV <- rbind(newV,center,a)
  }
}

# add new vertices to the sphere
if( nrow(newV) > 0 ){
  tess1 <- RefineSphericalTessellation( sphere$V, newV )
} else {
  tess1 <- sphere
}

# compute norming constant, possibly refine the tessellation, and get weights
# by computing the surface area of the contour
if( norm.const.method == "integrate" ) {
  f1 <- function( x ) { cfunc.eval( cfunc,  x )^cfunc$d }
  S1 <- aperm(tess1$S,c(2,1,3))
  a <- adaptIntegrateSphereTri( f1, S1, partitionInfo=TRUE, maxEvals=maxEvals, ... )
  if (a$returnCode != 0) {
    stop(paste("Error return from adaptIntegrateSphereTri: returnCode= ",a$returnCode," ",
      a$message, ".  Consider using norm.const.method='simplex.area'.", sep=""))
  }
  cfunc$norm.const <- 1.0/a$integral
  cfunc$tessellation.weights <- a$subsimplicesIntegral
  cfunc$functionEvaluations <- a$functionEvaluations
} else {
  if (norm.const.method=="simplex.area") {
    # just compute area of simplices
    n.tess1 <- dim(tess1$S)[3]
    cfunc$tessellation.weights <- double( n.tess1 )
    for (i in 1:n.tess1) {
      cfunc$tessellation.weights[i] <- SimplicialCubature::SimplexSurfaceArea(tess1$S[,,i])
    }
    cfunc$functionEvaluations <- 0    
    cfunc$norm.const <- 1/sum(cfunc$tessellation.weights)
    a <- list(subsimplices=tess1$S)
  } else {
    stop("Unknown value of norm.const.method: ",norm.const.method )
  }    
}

# construct a new mvmesh object using the partition info
tess2 <- mvmeshFromSimplices( aperm(a$subsimplices, c(2,1,3) ) )


# construct the  mvmesh for the contour by scaling points in tess2
V <- tess2$V
scale <- cfunc.eval( cfunc, t(V) )
for (i in 1:nrow(V)) {  
  V[i,] <- scale[i]*V[i,]
}
cfunc$tessellation <- mvmeshFromSVI( V, tess2$SVI, d-1 )

# add counts of the number of simplices at the three stages
cfunc$simplex.count <- c( ncol(sphere$SVI), ncol(tess1$SVI), ncol(tess2$SVI) )

class(cfunc) <- "gensphere.contour"
return(cfunc) }

#####################################################################
dgensphere <- function( x, gs.dist ){
# evaluate the generalized spherical function specified by gsfunc

stopifnot( class(gs.dist)=="gensphere.distribution" )
if( is.vector(x) ) { x <- matrix(x,ncol=1) }

v <- gs.vfunc.eval( gs.dist$cfunc, x )
p <- 1 - (gs.dist$cfunc)$d
z <- (gs.dist$cfunc)$norm.const * v^p * gs.dist$dradial( v )
j <- which(v==0)
if(length(j)>0) z[j] <- gs.dist$g0
return(z) }

#####################################################################
rgensphere <- function( n, gs.dist) {
# simulate n rarndom vectors from the generalized spherical
# two steps:
#  1. simulate from the contour/shell
#  2. scale each value from step 1 by a random amount given by the radial distribution

stopifnot( class(gs.dist)=="gensphere.distribution" )
x <- t( rmvmesh( n, (gs.dist$cfunc)$tessellation, (gs.dist$cfunc)$tessellation.weights ))
r <- gs.dist$rradial(n)
for (i in 1:n) {
  x[,i] <- r[i]*x[,i]
}
return(x)}

#####################################################################
gs.vfunc.eval <- function( cfunc, x ) {
#  evaluate v( x ) = |x| / c(x/|x|),  for x in {x[,i],i=1,...,ncol(x)},
#  where c(.)=cfunc.eval(cfunc,.) is a contour function described by cfunc

xnorm <- gs.lp.norm(x,p=2)
e1 <- double(nrow(x)); e1[1] <- 1.0  # first standard unit vector
u <- matrix(0.0,nrow=nrow(x),ncol=ncol(x))
for (i in 1:ncol(x)) {  # normalized each column of x[.,,]
  if (xnorm[i] > 0) { u[,i] <- x[,i]/xnorm[i] }
  else { u[,i] <- e1 }
}
cval <- cfunc.eval( cfunc, u )
return( xnorm/cval ) }

#####################################################################
gs.lp.norm <- function( x, p){
# compute the  l^p norm of the column vectors of x
y <- ( colSums( abs(x)^p ) )^(1.0/p)
return(y) }

#####################################################################
gs.gen.lp.norm <- function( x, p, A){
# compute the generalzed l^p norm of the column vectors of x, which
# is defined as gs.lp.norm( A %*% x, p), where
#    A is a (m x d) matrix
#    x is a (d x n) matrix, with columns being points in

# A is a (m x n) matrix, each column of A %*% x is an m vector
y <- gs.lp.norm( A %*% x, p ) 
return(y) }

#####################################################################
gs.elliptical <- function( x, B ) {
# compute sqrt(x^T B x) for each column of x
stopifnot( nrow(x)==nrow(B),nrow(B)==ncol(B) )
n <- ncol(x)
r <- double(n)
for (i in 1:n)  { 
  r[i] <- x[,i] %*% B %*% x[,i];
  #if(r[i] <0) cat(i,x[,i],"    ",r[i],"\n") 
}
return(sqrt(r))}

#####################################################################
gs.proj.normal <- function( x, mu, sigma ){
# computes a vector y, with
#         y[i] = |x[,i]| * exp( -z[i]/(2*sigma^2) )  or 0 depending on the 
# orientation of x[,i] and mu.  See details below.
#
#    x is a (d x n) matrix with columns x[,i] giving a point in R^d
#    mu is the "mean" vector, a point on the unit sphere in R^d
#    sigma > 0 is a scale, approximately the scale of the normal distribution.
#          using a small number approximates a point mass at mu;
#          using a large number approximates an indicator of the hemisphere
#              centered at mu.
#    z[i] = square of the distance between of mu and projection of
#          x[,i]/|x[,i]| onto the plane tangent to unit sphere at mu.
#    Notes
#    (1) Even when restricted to the unit sphere, this is not a normal
#        density because there is no norming constant, and it is not
#        easy to compute that constant due to the projection.
#    (2) This works in any dimension d > 1.
#    (3) The value is 0.0 if x[,i] is on the side of the plane
#        through the origin and perpendicular to mu that is opposite mu.
#    (4) 0 <= y[i] <= |x[,i]|.

k <- -1.0/(2.0*sigma*sigma)
mu.mu <- sum(mu^2)
m <- ncol(x)
y <- rep(0.0,m)
xnorm <- gs.lp.norm(x,p=2)
for (i in 1:m) {
  if (xnorm[i] > 0) {
    xx <- x[,i]/xnorm[i]  # unit vector in the direction of x[,i]
    a <- sum( xx*mu )
    if (a > 0) {  # return y[i]=0.0 if a <= 0
      xstar <- xx*mu.mu/a
      z <- sum( (xstar-mu)^2 )
      y[i] <- exp( k*z )
    }
  }
}
return(xnorm*y)}

#####################################################################
gs.cone <- function( x, mu, theta0 ){
# computes a vector y, with
#         y[i] = height of cone with peak 1 at center mu and height 0 r units
#                away from 0, where distance r from center is measured in radians
#                along the sphere.
#    x is a (d x m) matrix with columns x[,i] giving a unit vector in R^d
#    mu is the "mean" vector, a point on the unit sphere in R^d
#    0 < theta0 < pi/2 determines the location of the base of the cone:
#         all points whose angle = theta0
#    Notes
#    (1) This works in any dimension d > 1.
#    (2) The value is 0.0 if x[,i] is on the side of the plane
#        through the origin and perpendicular to mu that is opposite mu.
#    (3) 0 <= y[i] <= 1.

tmp <- min(pi/2,abs(theta0))
if (tmp != theta0) {
  warning(paste("theat0 changed to ",tmp," in function gs.cone"))
  theta0 <- tmp
} 
m <- ncol(x)
y <- double(m)
for (i in 1:m) {
  x.mu <- sum( x[,i]*mu )
  if( x.mu > 0 ) {
    if (x.mu > 1) { x.mu <- 1 }  # roundoff can make x.mu = 1+epsilon, 
                                # which causes acos(x.mu) to be NA
    theta <- acos(x.mu)
    ht <- 1.0 - theta/theta0
    if (ht > 0) { y[i] <- ht }
  }
}
return(y)}

#####################################################################
gs.pdf2d.plot <- function( gs.dist, xy.grid=seq(-10,10,.1) ){
# compute and plot the density surface z=dgensphere(x,y) on a grid

outer.f <- function( x, y, gs.dist ) { dgensphere( rbind(x,y), gs.dist) }
z <- outer( xy.grid, xy.grid, outer.f, gs.dist=gs.dist )

open3d()
rgl.surface(xy.grid,xy.grid,z,color="blue",back="lines")
aspect3d(1,1,1)
rgl.viewpoint( theta=15,phi=45 )
axes3d()

invisible(z) }

########################################################################
RefineSphericalTessellation <- function( V1, V2 ){
# define a new tessellation of the unit sphere in n-dimensions, using the
# tessellation obtained by merging the points in V1 and V2
# WARNING: it is assumed that one of V1 or V2 (or both) include
# the great circle of the sphere intersect {x[1]=1}.  This is true
# by construction if one of V1 or V2 is a UnitSphere object from
# the mvmesh package.

V <- rbind(V1,V2)
nV <- nrow(V)
n <- ncol(V1)

# separate into right and left hemispheres
# first treat the "right hemisphere"
which1 <- which( V[,1] >= 0 )

# project onto (n-1) dim. space
cur1 <- V[which1,2:n,drop=FALSE ]
tess1 <- delaunayn( cur1 )
# restore indices to original vertices
for (i in 1:nrow(tess1)) {
  for (j in 1:ncol(tess1)) {
    tess1[i,j] <- which1[tess1[i,j]]
  }
}

# repeat with the "left hemisphere", note this will include the
# great circle where x[1]=0 in both right and left hemisphere
which2 <- which( V[,1] <= 0 )
cur2<- V[which2,2:n,drop=FALSE ]
tess2 <- delaunayn( cur2 )   # shift to get position in V
for (i in 1:nrow(tess2)) {
  for (j in 1:ncol(tess2)) {
    tess2[i,j] <- which2[tess2[i,j]]
  }
}

# project both right and left parts onto sphere - don't need to
# do anything, just use original vertices in V and order in tess1 and tess2
SVI <- t( rbind( tess1, tess2 ) )

a <- mvmeshFromSVI( V, SVI, n-1 )
a$type <- "RefineSphericalTessellation"
return(a) }

########################################################################
NearbyPointsOnSphere <- function( x, epsilon ){
# add n=length(x) points approximately epsilon units away from x and equally
# spaced on the unit sphere.   New version allows epsilon to be a vector.
# it is assumed that x is a unit vector

stopifnot( is.vector(x), all(epsilon > 0), is.vector(x), length(x) > 1 )

n <- length(x)
# shift the unit simplex to be centered at (0,0,...,)
W1 <- diag( rep(1.0-1.0/n,n) )

# rotate to point in the direction x and shift to be centered at x
R <- RotateInNDimensions( rep(1/sqrt(n),n), x )

W2 <- matrix( 0.0, nrow=n*length(epsilon), ncol=n )
count <- 0
for (j in 1:length(epsilon)) {
  for (i in 1:n) {
    count <- count + 1
    W2[count, ] <- x + epsilon[j]*R %*% W1[i, ]
  }
}
# project onto unit sphere
for (i in 1:count) {
  norm <- sqrt( sum( W2[i, ]^2 ) )
  W2[i,] <- W2[i, ]/norm
}
return( W2 )}

########################################################################
RotateInNDimensions <- function( x, y ) {
# find an n x n matrix A that rotates vector x to the *direction* of the vector y
# if |x|=|y|, then A %*% x = y (up to numerical round-off)
# R implementation of the answer by Stephen Montgomery-Smith on the webpage
#    http://math.stackexchange.com/questions/598750/finding-the-rotation-matrix-in-n-dimensions

stopifnot( is.vector(x), is.vector(y), length(x)==length(y) )

norm.x <- sqrt(sum(x*x))
norm.y <- sqrt(sum(y*y))
u <- x/norm.x
v1 <- y - sum(u*y)*u
if ( sqrt(sum(v1^2) < 1.0e-14)) return( diag( rep(1.0,length(x) ) ) )
v <- v1/(sqrt(sum(v1*v1)))
cost <- sum(x*y)/(norm.x*norm.y)
sint <- sqrt( 1-cost^2)
Rt <- matrix( c(cost,sint,-sint,cost), 2, 2 )
uv <- cbind(u,v)
u <- matrix(u,nrow=length(u),ncol=1)
v <- matrix(v,nrow=length(v),ncol=1)
A <- diag( rep(1.0,length(u)) ) - (u %*% t(u)) - (v %*% t(v)) + uv %*% Rt %*% t(uv)
return(A)}

#######################################################################
#  printing and graphing functions for objects of class mvmesh
#######################################################################
print.gensphere.contour <- function( x, ... ) {
# quick method to print out contents of gensphere.contour x

if( class(x) != "gensphere.contour") warning( "x is not a gensphere.contour object" )

print(str(x),...) }

#######################################################################
plot.gensphere.contour <- function( x, multiplier=1, ... ) {
# quick method to plot a gensphere.contour x; the optional argument multiplier
# scales the contour by a factor of multiplier

if( class(x) != "gensphere.contour") warning( "x is not a gensphere.contour object" )

mesh <- x$tessellation
mesh$V <- multiplier*mesh$V
mesh$S <- multiplier*mesh$S
plot( mesh, ... ) }

#######################################################################
print.gensphere.distribution <- function( x, ... ) {
# quick method to print out contents of gensphere.distribution x

if( class(x) != "gensphere.distribution") warning( "x is not a gensphere.distribution object" )

print(str(x),...) }


