/*
Javier A. Porras Francisco
2019
*/

/**
Functions to assemble global matrices and vectors
from mesh coordinates or from local fragments.
*/
module femd.assembly;
@safe:
@nogc:
pure:
nothrow:

/**
Assembles the global mass matrix from 1D mesh `x`.

The mass matrix defines the linear equation
to obtain the orthogonal projection of any function
on the space of continuous piecewise linear functions
defined on mesh `x`.
The elements of the mass matrix are calculated
as the integrals across the mesh
(approximated by Simpson's formula)
of the multiplication in pairs of the "hat functions"
which form a convenient basis of this space.
See book §1.5.1.

Params:
M = (Output) assembled global mass matrix.
    All elements must be initialized to zero beforehand.
x = Mesh i.e. vector of node coordinates.
*/
void assembleMass1D(Matrix, Vector)
(ref Matrix M, const ref Vector x)
{
	size_t nnodes = x.length;
	for(size_t i = 1; i < nnodes; ++i) // nnodes-1 elements
	{
		size_t i1 = i - 1;
		auto
			h = x[i] - x[i1], // element's size.
			h3 = h/3,
			h6 = h/6;
		
		M[ i1 , i1 ] += h3;
		M[ i1 , i  ] += h6;
		M[ i  , i1 ] += h6;
		M[ i  , i  ] += h3;
	}
}
unittest ///
{
	import algebra.statics;

	// 3 nodes (2 elements) example:
	auto x = Vector!3( 0, 6, 18 );
	// This must be the result, easily calculated by hand:
	auto Mtest = Matrix!(3,3)(
		2, 1, 0,
		1, 6, 2,
		0, 2, 4);

	auto M = Mtest.zeros;
	assembleMass1D(M, x);
	assert(M == Mtest);

	// Translating the mesh has no effect on the mass matrix:
	x[] += 5;
	M = M.zeros;
	assembleMass1D(M, x);
	assert(M == Mtest);

	// Expanding the mesh multiplies the matrix:
	x[] *= 7;
	Mtest[] *= 7;
	M = M.zeros;
	assembleMass1D(M, x);
	assert(M == Mtest);

	// Translating again:
	x[] += 11;
	M = M.zeros;
	assembleMass1D(M, x);
	assert(M == Mtest);
}

/**
Assembles the global load vector from 1D mesh `x`.

The load vector of a function `f`
defines together with the mass matrix the linear equation
to obtain the orthogonal projection of `f`
on mesh `x`.
The elements of the load vector are calculated
as the integrals across the mesh
(approximated by the trapezoidal rule)
of `f` multiplied by each basis "hat function".
See book §1.5.2.

Params:
b = (Output) assembled global load vector.
    All elements must be initialized to zero beforehand.
x = Mesh i.e. vector of node coordinates.
f = Load function.
*/
void assembleLoad1D(Vector, scalar)
(ref Vector b, const ref Vector x, scalar delegate(scalar) pure nothrow @safe @nogc f)
{
	size_t nnodes = x.length;
	auto f1 = f(x[0]); // cache to avoid evaluating f() twice at every internal node
	for(size_t i = 1; i < nnodes; ++i) // nnodes-1 elements
	{
		size_t i1 = i - 1;
		auto
			x1 = x[i1],
			xi = x[i ],
			h2 = (xi - x1) / 2,
			fi = f(xi);
		
		b[i1] += h2 * f1;
		b[i ] += h2 * fi;

		f1 = fi;
	}
}
unittest ///
{
	real f(real) { return 42; } /* We choose a constant function because:
	a) it is already piecewise linear, so it should be equal to its projection; and
	b) it's a polynomial of degree 0, so its products with the hat functions
	are linear, and so the Trapezoidal approximation becomes exact. */
	
	import algebra.statics;
	// Mesh:
	auto x = Vector!3( -6, 0, 12 );
	// Mass matrix (by hand):
	auto M = Matrix!(3,3)(
		2, 1, 0,
		1, 6, 2,
		0, 2, 4);
	
	// Assemble load vector:
	auto b = x.zeros;
	assembleLoad1D(b, x, &f);

	// Check:
	assert( b == (M * Vector!3( f(x[0]), f(x[1]), f(x[2]) )).toVector );
}
