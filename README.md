# RSDMT
Random Strong Discrete Morse Theory

## Abstract
In topology, simplicial complexes can be reduced to smaller complexes of the same simple homotopy type through a series of elementary collapses. It is a common goal to try to find the smallest complex that is of the same simple homotopy type as another. We use the concept of strong collapses in which at least 2 faces are removed, to speed up the process of reducing the complexes. In our heuristic, we first check to see if we can run a strong collapse. If not, we attempt to perform a random elementary collapse if there are free pairs available. If there are no free pairs available, we then remove a random top face from the complex and try again. For the last two options, we keep track of how many top faces and elementary collapses were required in a vector. We obtain this strong discrete Morse vector by increasing the nth index in the vector each time an n-dimensional face is removed by a method other than a strong collapse. Since the heuristic relies on random selections, the vector may be different if run multiple times. The average vector obtained gives some indication of the complexity of the complex.

## Requirements
Below are the software requirements
  * Boost (that's it for now)
## Compile flags
  `-std=c++14`
## Etc.
[This](https://arxiv.org/pdf/1303.6422.pdf) paper is the inspiration for this project.
