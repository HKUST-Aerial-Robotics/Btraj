# Library design
We provide here a small overview of the software design choices.

### Criterions
- Code easy to use for the user.
- Avoid long and nested namespaces and variable names.
- Create a hierarchy that allows to reuse the code.
- Do not rely on weird C++ tricks.
- Do not rely on too many dependencies.
- Code flexibility and high performance.

### Design
- To reach high performance, `std::array` are mainly used.
- To avoid many dependencies, C++11 was chosen.
- To make the code easier to use, we only rely on Boost libraries. The other dependencies are already included in the code.
- To reuse code as much as possible, solvers are grouped depending on how they work. FMM-based solvers have many common procedures. The same as FM2-based solvers.
- Also, we tried to use the `information expert` pattern.
- In order to reach the flexibility and high performance, we chose to use a combination of static (templates) and dynamic  polymorphism. That is, all solvers inherit from the Solver class. However, Solver itself is templated depending the grid type.

This way, the compiler can carry out tons of optimizations by knowing in advance the number of dimensions on the grid. Actually, presumably the performance is the same as if the operations were hardcoded for any number of dimensions.

FMM-based solvers, following the `policies` design patter, have other parameter templates that change the behaviour. Concretely, the heap types.

The dynamic polymorphism allows to use all solvers under a common interface, as the examples included.

![Solvers hierarchy](solvers.png)