# Building
We tried to keep the CMake configuration as simple as possible. Therefore, there are 3 different things to build independently:

## main.cpp
Contains a demo of most of the main capabilities of the code. To compile and run:

	$ <go to fastmarching folder>
    $ mkdir build
    $ cd build
    $ cmake ..
    $ make
    $ ./fmm -map1 ../data/testimg.png -map2 ../data/map.png -vel ../data/velocities.png

`cmake` building type can be specified:

    $ cmake .. -DCMAKE_BUILD_TYPE=Release (Release, RelWithDebInf or Debug, Release by default)

## Examples
The examples folder contains three, more advanced examples. To compile them:

	$ <go to fastmarching folder>
    $ cd examples
    $ mkdir build
    $ cd build
    % cmake ..

As before, the cmake building type can be specified.

This will generate three binary files:

#### test_fmm.cpp
Contains a demo of running different FMM-based solvers:

	$ ./test_fmm

#### test_fm2.cpp
Contains a demo of running all the implemented FM2-based solvers:

	$ ./test_fm2 -map ../../data/map2.png

#### test_fmm_benchmark.cpp {#building_bmtest}
Contains a demo of a benchmark generated from code:

	$ ./test_fmm_benchmark

## Benchmarking application  {#building_bmapp}
Compiles a benchmark application which allows to automatically configure and run benchmarks by only creating a CFG file:

    $ <go to fastmarching folder>
    $ cd benchmark
    $ mkdir build
    $ cd build
    % cmake ..

Running a benchmark:

	$ ./fmm_benchmark ../benchmark.cfg

See the [detailed benchmark](markdown/benchmarking.md) documentation.
