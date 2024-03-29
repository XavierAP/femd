# femd
Under-construction self-study on
FEM (the Finite Element Method)
and the D programming language.

The original idea was to create a FEM library in D,
as a learning exercise in both.
In particular I take advantage of D's attractive features
in generic programming, constraints, compile-time evaluation, meta-programming, etc.
and language support for unit tests.

As a design choice, the [**femd**](/src/femd/) package requires linear algebra types,
but treats them in a fully generic, duck-typed way.
Therefore it's a template library with no dependencies
-- outside its unit tests.
To be able to write and run unit tests,
a minimum [**algebra**](/src/algebra/) package has been written,
and the unit tests, and only they, depend on it.

The [**matlab**](/matlab/) folder contains the Matlab implementations
copied from the book *The Finite Element Method: Theory, Implementation, and Applications* by Larson and Bengzon.
These implemenations are used to generate test benchmarks for the D re-implementations.
