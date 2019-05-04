/*
Javier A. Porras Francisco
2019
*/

/**
Simple Vectors and Matrices
of static (templated) dimensions
allocated on the stack.
*/
module algebra.statics;
@safe:
@nogc:
pure:
nothrow:

public
import algebra.scalar;

/**
Simple implementation of vectors of dimension fixed at compile time.
The default constructor returns an uninitialized (garbage-filled) object.
Params:
n  = Dimension, number of elements.
*/
struct Vector(size_t n)
{
	/// A vector is a tensor of order 1.
	static uint order() { return 1; }

	/// Returns the dimension.
	/// Returns 0 (instead of throwing) if dim is different than 1.
	static size_t length(uint dim = 1)
	{
		switch(dim)
		{
		case 1: return n;
		default: return 0;
		}
	}

	/// Indexed read.
	scalar opIndex(size_t pos) const
	{
		return mem[pos];
	}
	/// Indexed assignment.
	scalar opIndexAssign(scalar value, size_t pos)
	{
		return mem[pos] = value;
	}
	/// Indexed unary operators.
	scalar opIndexUnary(string op)(size_t pos)
	{
		return mixin(op~ "mem[pos]");
	}
	/// Indexed assignment operators.
	scalar opIndexOpAssign(string op)(scalar value, size_t pos)
	{
		return mixin("mem[pos]" ~op~ "= value");
	}
	

	/// Converts a Vector into a column Matrix.
	auto toColMatrix() const
	{
		return copyTo!(Matrix!(n,1));
	}
	/// Converts a Vector into a row Matrix.
	auto toRowMatrix() const
	{
		return copyTo!(Matrix!(1,n));
	}

	
	mixin common!(Vector!n); // Operators and methods common to Vectors and Matrices.

private:
	scalar[n] mem = void; /// Data storage.
	
}
private
unittest
{
	auto v = Vector!3();
	assert( v.length == 3 );
	assert( v.order == 1 );

	scalar a = 2;
	v[0] = a;
	assert( a == v[0] );

	assert( ++a == ++v[0] );
	assert(   a ==   v[0] );
	assert( --a == --v[0] );
	assert(   a ==   v[0] );

	v[0] += a;
	a += a;
	assert( a == v[0] );
	
	auto col = Matrix!(3,1)( 1, 2, 3 );
	auto row = Matrix!(1,3)( 1, 2, 3 );
	v = [ 1, 2, 3 ]; // assignment from array
	assert( v.toColMatrix == col );
	assert( v.toRowMatrix == row );
}


/**
Simple implementation of matrices of dimensions fixed at compile time.
The default constructor returns an uninitialized (garbage-filled) object.
Params:
nr = Number of rows.
nr = Number of columns.
*/
struct Matrix(size_t nr, size_t nc)
{
	/// A matrix is a tensor of order 2 (rows, columns).
	static uint order() { return 2; }

	/// Returns the dimension of the specified order (1=rows, 2=columns).
	/// Returns 0 (instead of throwing) if dim is different than 1 or 2.
	static size_t length(uint dim)
	{
		switch(dim)
		{
		case 1: return nr;
		case 2: return nc;
		default: return 0;
		}
	}

	/// Indexed read.
	scalar opIndex(size_t row, size_t col) const
	{
		return mem[row*nc + col];
	}
	/// Indexed assignment.
	scalar opIndexAssign(scalar value, size_t row, size_t col)
	{
		return mem[row*nc + col] = value;
	}
	/// Indexed unary operators.
	scalar opIndexUnary(string op)(size_t row, size_t col)
	{
		return mixin(op~ "mem[row*nc + col]");
	}
	/// Indexed assignment operators.
	scalar opIndexOpAssign(string op)(scalar value, size_t row, size_t col)
	{
		return mixin("mem[row*nc + col]" ~op~ "= value");
	}
	
	mixin common!(Matrix!(nr,nc)); // Operators and methods common to Vectors and Matrices.

	
	/// Converts a row or column Matrix into a Vector.
	auto toVector()() const
		if(nr == 1 || nc == 1)
	{
		return copyTo!(Vector!(mem.length));
	}

private:
	scalar[nc*nr] mem = void; /// Data storage.
	
}
private
unittest
{
	auto m = Matrix!(3,2)();
	assert( m.length(1) == 3 );
	assert( m.length(2) == 2 );
	assert( m.order == 2 );
	
	scalar a = 2;
	m[0,0] = a;
	assert( a == m[0,0] );

	assert( ++a == ++m[0,0] );
	assert(   a ==   m[0,0] );
	assert( --a == --m[0,0] );
	assert(   a ==   m[0,0] );

	m[0,0] += a;
	a += a;
	assert( a == m[0,0] );

	auto col = Matrix!(3,1)( 1, 2, 3 );
	auto row = Matrix!(1,3)( 1, 2, 3 );
	auto v = Vector!3( 1, 2, 3 );
	assert( v == col.toVector );
	assert( v == row.toVector );
}


/// Checks that two types (tensors of order 1 or 2) have the same size (e.g. 3d Vector and 3x1 or 1x3 Matrix).
template areSizesEqual(T1, T2)
{
	static if(is(T1 == T2))
		enum bool areSizesEqual = true;
	else
	{
		static if(T1.order == T2.order)
		{
			static assert(T1.length(1) != T2.length(1) || T2.length(2) != T2.length(2));
			enum bool areSizesEqual = false;
		}
		else static if(T1.order == 1)
		{
			static assert(T2.order == 2);
			enum bool areSizesEqual =
				(T2.length(1) == 1 && T2.length(2) == T1.length) ||
				(T2.length(2) == 1 && T2.length(1) == T1.length) ;
		}
		else static if(T1.order == 2)
		{
			static assert(T2.order == 1);
			enum bool areSizesEqual =
				(T1.length(1) == 1 && T1.length(2) == T2.length) ||
				(T1.length(2) == 1 && T1.length(1) == T2.length) ;
		}
		else static assert(false);
	}
}


private:

/// Operators and methods common to Vectors and Matrices.
mixin template common(T)
{
	/// Constructor with initialization.
	this(scalar[mem.length] elems...)
	{
		opAssign(elems);
	}
	/// Assignment operator.
	void opAssign(scalar[mem.length] elems)
	{
		mem = elems;
	}
	
	
	/// Factory method.
	static T zeros()
	{
		T ans;
		ans.mem[] = 0;
		return ans;
	}

	/// Sliced assignment operators.
	void opSliceOpAssign(string op)(scalar value)
	{
		mixin("mem[]" ~op~ "= value;");
	}

	/// Scalar post-multiplication.
	T opBinary(string op)(scalar x) const
		if(op == "*")
	{
		T ans;
		foreach(size_t k, elem; this.mem)
			ans.mem[k] = x * elem;
		
		return ans;
	}
	/// Scalar pre-multiplication.
	T opBinaryRight(string op)(scalar x) const
		if(op == "*")
	{
		return opBinary!op(x);
	}
	
	/// Addition and subtraction operators.
	T opBinary(string op, T2)(const auto ref T2 rhs) const
		if( (op == "+" || op == "-") && areSizesEqual!(T2, T) )
	{
		T ans = this;
		ans.opOpAssign!(op,T2)(rhs);
		return ans;
	}
	/// Assignment addition and subtraction operators.
	void opOpAssign(string op, T2)(const auto ref T2 rhs)
		if( (op == "+" || op == "-") && areSizesEqual!(T2, T) )
	{
		foreach(size_t k, ref elem; this.mem)
			mixin("elem " ~op~ "= rhs.mem[k];");
	}
	
	/// Matrix or Vector multiplication.
	auto opBinary(string op, T2)(const auto ref T2 rhs) const
		if(op == "*")
	{
		// 1. Compile-time check of dimensions and, in case T or T2 are Vectors, figure out whether to regard them as rows or columns.
		enum string errSizes = "Matrix multiplication dimension mismatch.";
		static if(T.order == 2)
		{
			enum size_t
				nr = T.length(1),
				n3 = T.length(2);
			
			static if(T2.order == 2)
			{
				static assert(n3 == T2.length(1), errSizes);
				enum size_t nc = T2.length(2);
			}
			else static if(T2.order == 1)
			{
				static if(n3 == 1)
				{
					enum size_t nc = T2.length;
				}
				else
				{
					static assert(n3 == T2.length, errSizes);
					enum size_t nc = 1;
				}
			}
			else static assert(false);
		}
		else static if(T.order == 1)
		{
			static if(T2.order == 2)
			{
				enum size_t
					n3 = T2.length(1),
					nc = T2.length(2);
				
				static if(n3 == 1)
				{
					enum size_t nr = T.length;
				}
				else
				{
					static assert(n3 == T.length, errSizes);
					enum size_t nr = 1;
				}
			}
			else static if(T2.order == 1)
			{
				// Scalar (dot) vector multiplication.
				enum size_t
					nr = 1,
					nc = 1,
					n3 = T.length;
				static assert(n3 == T2.length, errSizes);
			}
			else static assert(false);
		}
		else static assert(false); // order can only be 1 or 2.
		
		// 2. Algorithm. This naive one can't be outperformed except by multi-threading or in case of very huge matrices (which shouldn't be allocated statically as the types in this module).
		Matrix!(nr,nc) ans;
		foreach(size_t r; 0 .. nr)
		{
			foreach(size_t c; 0 .. nc)
			{
				size_t
					i = r*n3,
					j = c;
					
				scalar elem = 0;
				foreach(size_t k; 0 .. n3)
				{
					elem += this.mem[i+k] * rhs.mem[j];
					j += nc;
				}
				ans[r,c] = elem;
			}
		}
		return ans;
	}
	
	/// Infinity norm.
	scalar normInf() const
	{
		import std.math: abs;
		scalar ans = 0;
		foreach(elem; this.mem)
		{
			scalar ea = abs(elem);
			if(ea > ans) ans = ea;
		}
		assert( ans >= 0 );
		return ans;
	}

	private void copyFrom(T2)(const auto ref T2 rhs)
	{
		enum n = this.mem.length;
		static assert(n == rhs.mem.length);

		foreach(size_t k; 0 .. n)
			this.mem[k] = rhs.mem[k];
	}
	private T2 copyTo(T2)() const
	{
		T2 ans;
		ans.copyFrom(this);
		return ans;
	}
}


unittest
{
	auto v0 = Vector!2( 1, 2 );
	auto m0 = Matrix!(2,1)( 3, -1 );

	assert( v0 * m0 == Matrix!(1,1)( 1 ) );
	assert( m0 * v0 == Matrix!(2,2)( 3, 6, -1, -2 ) );

	assert( v0 == Vector!2.zeros + v0 );
	assert( m0 == Matrix!(2,1).zeros + m0 );

	scalar a = 8;
	auto v1 = v0;
	v1[] *= a;
	assert( v1 == a * v0 );
	assert( v1 == v0 * a );
	v1[] /= a;
	assert( v1 == v0 );

	auto m1 = Matrix!(1,2)( 0, -3 );
	assert( v0 + a * m0 - m1 == typeof(v0)(
		v0[0] + a * m0[0,0] - m1[0,0] ,
		v0[1] + a * m0[1,0] - m1[0,1] ) );

	v1 = v0;
	v1 += m0;
	v1 -= m0;
	assert( v1 == v0 );
}
