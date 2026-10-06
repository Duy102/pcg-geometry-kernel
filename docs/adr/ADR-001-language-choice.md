# ADR-001: Production language

**Decision:** C++20.

The current environment already provides GCC 14.2, Clang 17, CMake 3.31, Boost headers, GMP runtime support, and direct future interoperability with CGAL-style computational-geometry infrastructure. C++ also makes floating-point rounding control and Boost.Numeric.Interval readily available for the first certified vertical slice.

Rust remains attractive for memory safety, but it is not installed in the current execution environment and would add toolchain/bootstrap work before theorem validation. The project therefore chooses C++20 for the production kernel and compensates with a small API surface, value types, sanitizers in later phases, and strict tests.
