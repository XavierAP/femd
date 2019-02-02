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

See book 1.5 algorithm 2 (page 17).
The mass matrix contains the integrals of the products of the hat basis functions.
Simpson's formula is used for integration, see book 1.5.1.

Params:
M = (Output) assembled global mass matrix.
    All elements must be initialized to zero beforehand.
x = Vector of node coordinates.
*/
void assembleMass1D(Matrix, Vector)(ref Matrix M, const ref Vector x)
{
	size_t nnodes = x.length;
	for(size_t i = 1; i < nnodes; ++i) // nnodes-1 elements
	{
		size_t i1 = i - 1;
		auto h = x[i] - x[i1]; // element's size.
		
		auto h3 = h/3;
		auto h6 = h/6;
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
