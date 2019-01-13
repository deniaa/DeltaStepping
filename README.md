# DeltaStepping
Code contains simple Dijkstra algorithms (single thread and parallel mode) for compare results with delta-stepping. Single thread Dijkstra execution can take alot of time on huge graphs. Sorry for no option about this (temporary) misunderstanding.

Single thread simple delta-stepping algorithm really work correctly and faster even than the parallel Dijkstra algorithm on graph with small maximum vertex degree.

Improved delta stepping algorithm probably contains a bug... because of result is not equal with results of other algorithms.

## Compiling
Required /openmp key in compilation for working really in parallel mod in some algorithms.
For windows native compiler: ```cl delta.cpp /openmp```

## Execution
Program have two required argumets: vertex count and maximum degree of vertex.
For windows: ```delta.exe <vertex_count> <maximum_degree>``` (for example ```delta.exe 1000 10```)

Real difference in execution time between the algorithms noticed on big graphs with small maximum vertex degree 
For example on 100000 vertexes and 50 maximum degree on my machine:
```
>delta.exe 100000 50
Graph generating time: 0.185342
6 threads Dijkstra calculating time: 1.723321
Single thread Dijkstra calculating time: 8.420936
Delta-stepping calculation time: 0.327579
Delta-stepping improved parallel (6 threads) calculation time: <Infinite loop>
```
