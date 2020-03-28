
# Exploratory/play code to figure out distribution of site sizes.
#
# This checks that the mixture of two uniforms has the correct properties in
# terms of mean and variance.
#
# Used to modify the random data generator code in blkvar package.

calc.rat = function( N ) {
  
  R = 100000
  X = 50
  N = N * X
  p = (N-X)/N
  A = runif( R*p, 0, X )
  B = runif( R*(1-p), X, N )
  Y = c( A, B )
  
  mean( Y )
  sd( Y )
  
  data.frame( EY = mean(Y), sd.Y = sd( Y ), rat = sd( Y ) / mean( Y ) )
}


calc.rat( 2 )
calc.rat( 3 )
calc.rat( 4 )

mults = seq( 1.1,10, length.out=200 )
rats = map_df( mults, calc.rat )
head( rats )
rats$mult = mults
qplot( mults, rats$rat )

qplot( mults, rats^2 )

qplot( Y )

lm( rats^2 ~ mults )

rat2 = rats$rat^2
lm( mults ~ rat2 )


-1.867 + 6.365 * 1/3

calc.rat( 2 )

calc.rat( 2 )^2


1 + 3 * 0.6165607


block.distn = function( J, n.bar, size.ratio ) {
  N = 1 + 3 * size.ratio
  p = (N-1)/N
  small = rbinom( J, 1, p )
  Y = runif( J )
  Y = n.bar * ifelse( small, Y, Y*(N-1) + 1 )
  
  Y
}

block.distn( 20, 50, 1/3 )


Y = block.distn( 100000, 50, 1/3 )
mean( Y )
sd( Y )
var(Y) / mean( Y )^2
