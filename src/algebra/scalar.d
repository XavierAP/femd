/*
Javier A. Porras Francisco
2019
*/

/**
Defines the scalar type used throughout the algebra package.
*/
module algebra.scalar;


version(PRECISION64)
	alias scalar = double;
else
	alias scalar = real;
