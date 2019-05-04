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


/**
Assembles the stiffness matrix from 1D mesh `x`
and Robin boundary conditions at both ends of the mesh.

The stiffness matrix K+R defines
together with the load vector b+r
the finite-element form (K+R)u=b+r of the problem
	-( a(x)·u'(x) )' = f(x)
with the Robin boundary conditions
	a(x0)·u'(x0) = k0·( u(x0) - g0 )
	a(xN)·u'(xN) = kN·( u(xN) - gN )

The terms R and r arise from the most general Robin boundary conditions.
Any Dirichlet or Neumann condition can be expressed
as a particular case of a Robin condition.

See book §2.4.1.

Params:
KR = (Output) assembled global stiffness matrix.
     All elements must be initialized to zero beforehand.
x  = Mesh i.e. vector of node coordinates.
a  = Stiffness function (see equations in Description).
     Must be strictly positive everywhere on `x`.
k0 = Boundary condition coefficient.
     Must be >= 0.
kN = Boundary condition coefficient.
     Must be >= 0.
*/
void assembleStiffness1D(Matrix, Vector, scalar)
(ref Matrix KR, const ref Vector x,
scalar delegate(scalar) pure nothrow @safe @nogc a,
scalar k0, scalar kN)
{
	size_t n = x.length - 1;
	for(size_t i = 1; i <= n; ++i)
	{
		size_t i1 = i - 1;
		auto amh = a((x[i] + x[i1])/2) / (x[i] - x[i1]);
		KR[ i1, i1 ] += amh;
		KR[ i1, i  ] -= amh;
		KR[ i , i1 ] -= amh;
		KR[ i , i  ] += amh;
	}
	KR[0,0] += k0;
	KR[n,n] += kN;
}
unittest ///
{
	/* Steady-state temperature in a rod, from book §2.4.2. */
	
	real conduct(real x) { return 0.5 - 0.06 * x; } // conductivity
	real // boundary conditions:
		k0 = 1e6,
		kN = 0;
	
	import algebra.statics;
	enum size_t n = 6;
	Vector!n x = [ 2, 2.2, 2.4, 2.6, 2.8, 3 ];
	auto K = Matrix!(n,n).zeros;
	
	assembleStiffness1D(K, x, &conduct, k0, kN);
	
	assert( 1e-9 > (K - Matrix!(n,n)( k0 +
		 1.87, -1.87,  0,     0,     0,     0,
		-1.87,  3.68, -1.81,  0,     0,     0,
		 0,    -1.81,  3.56, -1.75,  0,     0,
		 0,     0,    -1.75,  3.44, -1.69,  0,
		 0,     0,     0,    -1.69,  3.32, -1.63,
		 0,     0,     0,     0,    -1.63,  1.63
		+ kN )).normInf );
}

/**
Assembles the load vector from 1D mesh `x`
and Robin boundary conditions at both ends of the mesh.

See `assembleStiffness1D` for math details.

Params:
br = (Output) assembled global load vector.
     All elements must be initialized to zero beforehand.
x  = Mesh i.e. vector of node coordinates.
f  = Load function.
k0 = Boundary condition coefficient.
     Must be >= 0.
g0 = Boundary condition independent term.
kN = Boundary condition coefficient.
     Must be >= 0.
gN = Boundary condition independent term.
*/
void assembleLoad1D(Vector, scalar)
(ref Vector br, const ref Vector x,
scalar delegate(scalar) pure nothrow @safe @nogc f,
scalar k0, scalar g0, scalar kN, scalar gN)
{
	assembleLoad1D(br, x, f);
	size_t n = x.length - 1;
	br[0] += k0 * g0;
	br[n] += kN * gN;
}
unittest ///
{
	/* Steady-state temperature in a rod, from book §2.4.2. */
	
	import std.math: pow;
	real f(real x) { return 0.03 * (x - 6).pow(4); } // heat generation
	real // boundary conditions:
		k0 = 1e6,
		g0 = -1,
		kN = 0,
		gN = 0;
	
	import algebra.statics;
	enum size_t n = 6;
	Vector!n
		x = [ 2, 2.2, 2.4, 2.6, 2.8, 3 ],
		b = x.zeros;
	
	assembleLoad1D(b, x, &f, k0, g0, kN, gN);
	
	assert( 1e-9 > (b - Vector!n( k0 * g0 +
		0.768, 1.2510816, 1.0077696, 0.8018016, 0.6291456, 0.243
		+ kN * gN )).normInf );
}
