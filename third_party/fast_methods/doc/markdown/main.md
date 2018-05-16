These pages provide basic documentation about the Fast Marching code packages I have created.

- [Software design](markdown/design.md)
- [Detailed building](markdown/building.md)
- [Benchmarking](markdown/benchmarking.md)
- [nDGridMap implementation](ndgridmap.pdf): Probably, the most complex part of the code. It encapsulates the complexity of the n-dimensional generalization, so that solvers are implemented almost independently of the number of dimensions.

Please, carefully read the [README](https://github.com/jvgomez/fastmarching/blob/master/README.md/).

\note Some of the classes are just wrappers that can be used with lower level classes. For instance, declaring a `SFMM<nDGridMap<FMCell,2>>` is the same as using `FMM<nDGridMap<FMCell,2>, FMPriorityQueue<FMCell>>`. The same happens with the *Star solvers.
 
\note Due to the nDimensional implementation, the Eikonal solver code for a given index is probably not very optimized, since it computes from 1D to nD until a better solution is not possible. For 2D and 3D (and probably nD) this is not necessary, but my tests while implementing FSM obligated me to do it that way in order to have consistent results.
