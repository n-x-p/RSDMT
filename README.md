# RSDMT
Random Strong Discrete Morse Theory

## Abstract
In topology, simplicial complexes can be reduced to smaller complexes of the same simple homotopy type through a series of elementary collapses. It is a common goal to try to find the smallest complex that is of the same simple homotopy type as another. We use the concept of strong collapses in which at least 2 faces are removed, to speed up the process of reducing the complexes. In our heuristic, we first check to see if we can run a strong collapse. If not, we attempt to perform a random elementary collapse if there are free pairs available. If there are no free pairs available, we then remove a random top face from the complex and try again. For the last two options, we keep track of how many top faces and elementary collapses were required in a vector. We obtain this strong discrete Morse vector by increasing the nth index in the vector each time an n-dimensional face is removed by a method other than a strong collapse. Since the heuristic relies on random selections, the vector may be different if run multiple times. The average vector obtained gives some indication of the complexity of the complex.

## Requirements
Below are the software requirements
  * Boost >= 1.58.0
  * (that's it for now)
## Compile flags
  `-std=c++14`
## Possibly Confusing Stuff
The simplicial complex is stored in a boolean sparse matrix that represents a the adjacency matrix of a [Hasse diagram](https://en.wikipedia.org/wiki/Hasse_diagram) to make the program lighter on system memory.
## Inspiration
[This](https://arxiv.org/pdf/1303.6422.pdf) paper, **Random Discrete Morse Theory and a New Library of Triangulations**, by Bruno Benedetti and Frank H. Lutz is the inspiration for this project.
