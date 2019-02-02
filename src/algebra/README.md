# algebra package

## [**statics**](/src/algebra/statics.d) module
This module provides a simple implementation of Matrices and Vectors
intended for lightweight use
(e.g. not very large or sparse; static sizes).
The two type templates
Matrix and Vector
have the following features:
+ Matrices and Vectors of different sizes are different types (structs)
with no inheritance relationship;
however their operators are compatible when sizes are;
this is checked at compile time.
+ Operations:
indexing,
addition,
subtraction,
multiplication,
scalar pre-multiplication
and scalar assignment operators.
+ Stack allocation,
no garbage collector.
+ Matrices and Vectors may be constructed in one of three ways:
	+ Specifying every element.
	+ Static factory methods, e.g. `zeros()`.
	+ Default constructors: uninitialized, start filled with garbage.
