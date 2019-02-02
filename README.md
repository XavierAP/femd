# femd
Under-construction self-study on
FEM (the Finite Element Method)
and the D programming language.

The original idea was to create a FEM library in D,
as a learning exercise in both.
In particular I take advantage of D's attractive features
in generic programming, constraints, compile-time evaluation, meta-programming, etc.

As a design choice, the [**femd package**](/src/femd/) requires linear algebra types,
but treats them in a fully generic, duck-typed way.
Therefore it's a template library with no dependencies
-- outside its unit tests.
To be able to write and run unit tests,
a minimum [**algebra package**](/src/algebra/) has been written,
and the unit tests, and only they, depend on it.
