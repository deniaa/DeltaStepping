# DeltaStepping

## Compiling
Required /openmp key in compilation for working really in parallel mod in some algorithms.
For windows native compiler: ```cl delta.cpp /openmp```

## Execution
Program have two required argumets: vertex count and maximum degree of vertex.
For windows: ```delta.exe <vertex_count> <maximum_degree>``` (for example ```delta.exe 1000 10```)
