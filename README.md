Molekularis is a puzzle type invented by Alexander Thomas. Puzzles can
be found on his [website](http://molekularis.free.fr/). This project
provides a generator/solver for these types of puzzles as well as gui
for humans to solve them. You can try them online
[here](https://smanecke.github.io/molekularis.html).

 # Rules

Atoms of type H, O, N, C are layed out on a hexagonal grid. To solve
the puzzle, add a number of connections (from zero to three) between
neighboring atoms, such that

1. The number of connections to an atom of type H, O, N or C are one, two, three, and four, respectively.
2. Every two atoms are connected by a path formed by the connections.

# Build instructions

To build, run `make`. The `Makefile` is simple enough to be edited by hand, if needed.
The default compiler is `clang`, as it is the only major C-compiler which targets
`wasm`.

The `gui` program uses raylib. Uncomment the corresponding line in the Makefile to build it.

The Makefile builds into a seperate directory `build`.

# Usage

## Generator

The *generator* can be run by invoking `./build/molecularis size [#H #O #N #C]`
where size is the name of a file containing a template describing how
the puzzle should look like. We provide the templates 'tiny', 'small',
'medium', 'large', 'gigantic'. The integers #H, #O, #N, #C indicate
the probabilities that the respective atoms are chosen. The program 
generates a puzzle with an *unique solution* of the respective 
size and writes it into 'puzzle.txt'.

For example,

```
./build/molecularis small 0 1 10 10
```

generates a (small) puzzle where __H-atoms__ are never chosen and
__O-atoms__ 10 times less likely as __N-__ or __C-atoms__.
Superficially, a puzzle gets harder if it has very few __H-atoms__ and
__O-atoms__. Conversely the generator takes longer on these inputs, as
the probability of a unique solution for a random puzzle is very small.

## Gui

The *gui* solver can be run by invoking `./build/gui`. It reads to
puzzle in `puzzle.txt` and displays it for a human to solve it. The
left and right mouse button adds or removes a bond. Middle-clicking a
particular bond denotes that this count is finishes. Pressing __z__ on
the keyboard reverts the last step.

The gui does not check if the puzzle is correctly solved.

## Web
The web version can be run by pointing a web server to the `html` directory and pointing a browser to it. For example

```
python3 -m http.server --directory html
```

The web version uses the same controls as the gui. Progress is
automatically saved locally, so the puzzle does not need to be solved
in one sitting.

# How it works

The generator is very simple: try a random assignment and check if it
has a unique solution. If not, that is, if it either has no solutions
or at least two, do a small change and try again. So the interesting
part is the solver. The problem is translated into a SAT program using
standard techniques (see, for example, Donald Knuth's chapter on SAT
solving in `The art of programming`). But, this is only done for the
first of the two requirements listed in the rules. The connectivity
requirement is added in afterwards, by adding further constraints if
the solution is not connected.

This is done as follows: Assume that the solution has two connected
components. Then there are some some bonds which would connect them,
but which were not provided in the solution. In graph theory parlance,
these would be called *cut edges*. So one of these cut edges needs to
have at least one bond, but this condition is trivial to formulate as
a SAT clause. Whenever we find that the solution is not connected, we
add a clause corresponding to the cut edges and try again until we
either have a correct solution or found that the puzzle is not
solvable.

We offload the solving part to picosat, which is provided as an
`.c/.h` pair in the sources. 