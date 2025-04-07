
# Centrality Minimization in Graphs

This program computes and minimizes **betweenness** and **stress centrality** in
undirected graphs by evaluating the impact of adding missing edges.

## Compilation

To compile the code, run:

```bash
g++ -fopenmp final.cpp -o final
```

Make sure you have a C++ compiler with OpenMP support.

## Usage

```bash
./final <input_file>
```

If no input file is provided, the program will prompt you to enter one manually.

## Input Format

The input graph should be in a suitable format expected by `readGraph(filename)`.

## Output

The program outputs:

- Initial and minimized centrality values.
- Edges responsible for minimum betweenness and stress.
- Runtime performance.
