#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <set>
#define ARRSIZE(x)  (sizeof(x) / sizeof((x)[0])

double delta = 0.1;
struct edge
{
    double weight;
    int vertex;
};

struct weghtedGraph
{
    std::vector<struct edge>* graph;
    int* sizes;    
};

struct request
{
	int vertex;
	double weight;
};

struct vertex
{
	int u;
	int v;
	double w;
	
	bool operator < (const vertex& v) const
    {   
        if (this->u > v.u)
            return true;
        if (this->u < v.u)
            return false;
		
        if (this->v > v.v)
            return true;
        if (this->v < v.v)
            return false;
		
        if (this->w > v.w)
            return true;
		
        return false;
    }
	
    bool operator == (const vertex& v) const
    {
        if (this->u== v.u && this->v == v.v && this->w == v.w)
            return true;
        return false;
    }
};

int hashSize;
std::vector<struct vertex>* hash;

int GetHashIndex(int u, int v)
{
	return (u % hashSize) * hashSize + v % hashSize;
}

void PutToHash(int u, int v, double w)
{
	struct vertex vtx;
	vtx.u = u;
	vtx.v = v;
	vtx.w = w;
	hash[GetHashIndex(u, v)].push_back(vtx);
}

double GetFromHash(int u, int v)
{
	std::vector<struct vertex> bucket = hash[GetHashIndex(u, v)];
	
	for (int i = 0; i < bucket.size(); i++)
	{
		if (bucket[i].u == u && bucket[i].v == v)
			return bucket[i].w;
	}
	
	return INT_MAX;
}

void InitializeHash(int size)
{
	hashSize = size;
	hash = new std::vector<struct vertex>[size * size];
}

void FindShortestPathes(int vertexCount, int startVertex, struct weghtedGraph wGraph, int* parent, double* pathWeight)
{
    int* visited = new int[vertexCount];
    for (int i = 0; i < vertexCount; i++)
    {
        visited[i] = 0;
        pathWeight[i] = INT_MAX;
        parent[i] = -1;
    }
    pathWeight[startVertex] = 0;
    
    int* sizes = wGraph.sizes;
    std::vector<struct edge>* graph = wGraph.graph;
    
    for (int i = 0; i < vertexCount; i++)
    {
        int next = -1;
        int j;
        for (j = 0; j < vertexCount; j++)
        {
            if (!(visited[j] == 1) && (next == -1 || pathWeight[j] < pathWeight[next]))
                next = j;
        }
           
        if (pathWeight[next] == INT_MAX)
            break;
        visited[next] = 1;
            
        int count = sizes[next];
        for (j = 0; j < count; j++)
        {
            int to = graph[next][j].vertex;
            double len = graph[next][j].weight;
            if (pathWeight[next] + len < pathWeight[to]) {
                pathWeight[to] = pathWeight[next] + len;
                parent[to] = next;
            }
        }
    }
    
    free(visited);
}

void FindShortestPathesParallel(int vertexCount, int startVertex, struct weghtedGraph wGraph, int* parent, double* pathWeight)
{
    int* visited = new int[vertexCount];
    for (int i = 0; i < vertexCount; i++)
    {
        visited[i] = 0;
        pathWeight[i] = INT_MAX;
        parent[i] = -1;
    }
    pathWeight[startVertex] = 0;
    
    int* sizes = wGraph.sizes;
    std::vector<struct edge>* graph = wGraph.graph;
    
    //private variables
    int threadId;
    int threadsCount;
    int firstIndex;
    int lastIndex;
    
    int firstIndexForUpdate;
    int lastIndexForUpdate;
    int countForUpdate;
    
    
    //shared variables
    double minimumDistance;
    int nextVertex;
    
    # pragma omp parallel private (threadId, threadsCount, firstIndex, lastIndex, firstIndexForUpdate, lastIndexForUpdate)  shared (visited, pathWeight, graph, sizes, minimumDistance, nextVertex, countForUpdate)
    {
        threadId = omp_get_thread_num();
        threadsCount = omp_get_num_threads(); 
        firstIndex = (threadId * vertexCount) / threadsCount;
        lastIndex = ((threadId + 1) * vertexCount) / threadsCount - 1;
        
        for (int i = 0; i < vertexCount; i++)
        {
            # pragma omp single 
            {
                minimumDistance = INT_MAX;
                nextVertex = -1;
            }
            
            int localNextVertex = -1;
            double localMinimumDistance = INT_MAX;
            int j;
            for (j = firstIndex; j <= lastIndex; j++)
            {
                if (!(visited[j] == 1) && (localNextVertex == -1 || pathWeight[j] < localMinimumDistance))
                {
                    localNextVertex = j;
                    localMinimumDistance = pathWeight[j];
                }
            }
            
            # pragma omp critical
            {
                if (localMinimumDistance < minimumDistance)  
                {
                    minimumDistance = localMinimumDistance;
                    nextVertex = localNextVertex;
                }
            }
            
            # pragma omp barrier
            
            if (nextVertex == -1)
                break;
            
            # pragma omp single 
            {
                visited[nextVertex] = 1;                
                countForUpdate = sizes[nextVertex];                
            }
            
            # pragma omp barrier
            
            
            firstIndexForUpdate = (threadId * countForUpdate) / threadsCount;
            lastIndexForUpdate = ((threadId + 1) * countForUpdate) / threadsCount - 1;
                
          
            for (j = firstIndexForUpdate; j <= lastIndexForUpdate; j++)
            {
                struct edge edg = graph[nextVertex][j];
                int to = edg.vertex;
                double newLength = pathWeight[nextVertex] + edg.weight;
                if (newLength < pathWeight[to]) 
                {
                    pathWeight[to] = newLength;
                    parent[to] = nextVertex;
                }
            }
            
            #pragma omp barrier
        }
    }
    
    free(visited);
}

enum Kind { light, heavy };

bool wayToSort(struct edge const& i, struct edge const& j) 
{ 
	return i.weight < j.weight; 
}

struct weghtedGraph GenerateGraph(int n, int maxDegree)
{
    struct weghtedGraph g;
    g.graph = new std::vector<struct edge>[n];
    int* sizes = new int[n];
    
    for (int i = 0; i < n; i++)
        sizes[i] = 0;
    
    srand ( time ( NULL));
	#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        int vertexDegree = (maxDegree) * ( (double)rand() / (double)RAND_MAX );
		if (vertexDegree <= 0)
			vertexDegree = 1;
		std::vector<struct edge> neighbours;
		g.graph[i] = neighbours;
		//printf("%d'th vertex will have %d neighbours:\n", i, vertexDegree);

        for (int j = 0; j < vertexDegree; j++)
        {
            double weight = (double)rand()/RAND_MAX;
			int next;
			while (1)
			{
				next = (n) * ( (double)rand() / (double)RAND_MAX );
				int correct = 1;
				if (next == i)
					correct = 0;
				if (correct == 1)
				{
					for (int k = 0; k < j; k++)
					{
						if (g.graph[i][k].vertex == next)
						{
							correct = 0;
							break;
						}
					}
				}
				
				if (correct)
				{
					//printf("\t This is %d with weight %f\n", next, weight);
					struct edge e;
					e.vertex = next;
					e.weight = weight;
					g.graph[i].push_back(e);
					sizes[i] = sizes[i] + 1;
					break;
				}
			}
			
			//printf("%d\n", g.graph[i].size());
        }
		
		std::sort(g.graph[i].begin(), g.graph[i].end(), wayToSort);
    }
    
    g.sizes = sizes;
	
    return g;
}

void relax(std::vector<std::set<int>> &buckets, double* tent, struct request req)
{
	//std::cout << "\nNext relax - weight: " << req.weight << " Vertex: " << req.vertex << " tent: " << tent[req.vertex] << std::endl;

	if (req.weight < tent[req.vertex])
	{
		//printf("In if\n");
		
		if (tent[req.vertex] != INT_MAX)
		{
			int index = ((int)(tent[req.vertex] / delta)) % buckets.size();
			buckets[index].erase(req.vertex);			
			//printf("Erasing from %d bucket %d vertex\n", index, req.vertex);
		}
		
		int index = ((int)(req.weight / delta)) % buckets.size();
		buckets[index].insert(req.vertex);
		//printf("Inserting: bucket[%d] size = %zd\n", index, buckets[index].size());
		
		tent[req.vertex] = req.weight;
		//printf("Updating tent\n");
	}
}

void findRequests(struct weghtedGraph weightedGraph, double* tent, std::set<int> &bucket, Kind kind, std::vector<struct request> &requests)
{
	std::set<int>::iterator it;
	for (it = bucket.begin(); it != bucket.end(); ++it)
	{
		int vertex = *it;
		std::vector<struct edge> neighbours = weightedGraph.graph[vertex];
		int neighboursCount = weightedGraph.sizes[vertex];
		if (kind == Kind::light)
		{
			for (int i = 0; i < neighboursCount; i++)
			{
				struct edge neighbour = neighbours[i];
				if (neighbour.weight <= delta)
				{
					struct request req;
					req.vertex = neighbour.vertex;
					req.weight = tent[vertex] + neighbour.weight;
					requests.push_back(req);
				}
				else
				{
					break;
				}
			}
		}
		else
		{
			for (int i = neighboursCount - 1; i >= 0; i--)
			{
				struct edge neighbour = neighbours[i];
				if (neighbour.weight > delta)
				{
					struct request req;
					req.vertex = neighbour.vertex;
					req.weight = tent[vertex] + neighbour.weight;
					requests.push_back(req);
				}
				else
				{
					break;
				}
			}
		}
	}
}

void findRequestsParallel(struct weghtedGraph weightedGraph, double* tent, std::set<int> &bucket, Kind kind, std::vector<struct request> &requests)
{
	int bucketSize = bucket.size();
	int* tmp = new int[bucketSize];
	std::set<int>::iterator it;
	int ind = 0;
	for (it = bucket.begin(); it != bucket.end(); ++it)
	{
		tmp[ind] = *it;
		//printf("%d in tmp at index %d\n", *it, ind);
		ind++;
	}
	
	int j = 0;
	# pragma omp parallel for shared (tmp)
	for (j = 0; j < bucketSize; j++)
	{	
		std::vector<struct request> localRequests;
		
		int vertex = tmp[j];
		std::vector<struct edge> neighbours = weightedGraph.graph[vertex];
		int neighboursCount = weightedGraph.sizes[vertex];
		//printf("%d (parallel) neighbours of %d at %d index\n", neighboursCount, vertex, j);
		if (kind == Kind::light)
		{
			for (int i = 0; i < neighboursCount; i++)
			{
				struct edge neighbour = neighbours[i];
				if (neighbour.weight <= delta)
				{
					struct request req;
					req.vertex = neighbour.vertex;
					req.weight = tent[vertex] + neighbour.weight;
					localRequests.push_back(req);
				}
				else
				{
					break;
				}
			}
		}
		else
		{
			for (int i = neighboursCount - 1; i >= 0; i--)
			{
				struct edge neighbour = neighbours[i];
				if (neighbour.weight > delta)
				{
					struct request req;
					req.vertex = neighbour.vertex;
					req.weight = tent[vertex] + neighbour.weight;
					localRequests.push_back(req);
				}
				else
				{
					break;
				}
			}
		}
		
			
		#pragma omp critical
		{
			for (int i = 0; i < localRequests.size(); i++)
			{
				requests.push_back(localRequests[i]);
			}
		}
	}
}

void FindShortestPathesDeltaStepping(int vertexCount, int startVertex, struct weghtedGraph wGraph, double* tent)
{
	for (int i = 0; i < vertexCount; i++)
		tent[i] = INT_MAX;
	
	int bucketsCount = (double)1 / delta + 1;
	std::vector<std::set<int>> buckets (bucketsCount);
	
	
	struct request initialRequest;
	initialRequest.vertex = 0;
	initialRequest.weight = 0;
	relax(buckets, tent, initialRequest);
	int cyclingIndex = 0;
	while (true)
	{
		int minIndex = -1;
		for (int i = 0; i < bucketsCount; i++)
		{
			int indexForCheck = i + cyclingIndex;
			if (indexForCheck >= bucketsCount)
			{
				indexForCheck -= bucketsCount;
			}
			if (buckets[indexForCheck].size() > 0)
			{
				minIndex = i;
				break;
			}
		}
		cyclingIndex++;
		if (cyclingIndex >= bucketsCount)
			cyclingIndex = 0;
		
		if (minIndex == -1)
		{
			break;
		}		

		std::set<int> currentSet;
		
		while(!buckets[minIndex].empty())
		{
			std::vector<struct request> requests;
			findRequests(wGraph, tent, buckets[minIndex], Kind::light, requests);
			std::set<int>::iterator it;
			for (it = buckets[minIndex].begin(); it != buckets[minIndex].end(); ++it)
			{
				int vertex = *it;
				currentSet.insert(vertex);
			}
			buckets[minIndex].clear();
			
			for (int i = 0; i < requests.size(); i++)
			{
				relax(buckets, tent, requests[i]);
			}
			requests.clear();
		}
		
		
		std::vector<struct request> heavyRequests;
		findRequests(wGraph, tent, currentSet, Kind::heavy, heavyRequests);
		for (int i = 0; i < heavyRequests.size(); i++)
		{
			relax(buckets, tent, heavyRequests[i]);
		}
		heavyRequests.clear();
	}
}

void FindShortestPathesDeltaSteppingParallel(int vertexCount, int startVertex, struct weghtedGraph wGraph, double* tent)
{
	for (int i = 0; i < vertexCount; i++)
		tent[i] = INT_MAX;
	
	int bucketsCount = (double)1 / delta + 1;
	std::vector<std::set<int>> buckets (bucketsCount);
	
	
	struct request initialRequest;
	initialRequest.vertex = 0;
	initialRequest.weight = 0;
	relax(buckets, tent, initialRequest);
	int cyclingIndex = 0;
	while (true)
	{
		int minIndex = -1;
		for (int i = 0; i < bucketsCount; i++)
		{
			int indexForCheck = i + cyclingIndex;
			if (indexForCheck >= bucketsCount)
			{
				indexForCheck -= bucketsCount;
			}
			if (buckets[indexForCheck].size() > 0)
			{
				minIndex = i;
				break;
			}
		}
		cyclingIndex++;
		if (cyclingIndex >= bucketsCount)
			cyclingIndex = 0;
		
		if (minIndex == -1)
		{
			break;
		}		

		std::set<int> currentSet;
		
		while(!buckets[minIndex].empty())
		{
			std::vector<struct request> requests;
			findRequestsParallel(wGraph, tent, buckets[minIndex], Kind::light, requests);
			
			std::set<int>::iterator it;
			for (it = buckets[minIndex].begin(); it != buckets[minIndex].end(); ++it)
			{
				int vertex = *it;
				currentSet.insert(vertex);
			}
			buckets[minIndex].clear();
			
			for (int i = 0; i < requests.size(); i++)
			{
				relax(buckets, tent, requests[i]);
			}
			requests.clear();
		}
		
		
		std::vector<struct request> heavyRequests;
		findRequestsParallel(wGraph, tent, currentSet, Kind::heavy, heavyRequests);
		
		for (int i = 0; i < heavyRequests.size(); i++)
		{
			relax(buckets, tent, heavyRequests[i]);
		}
		heavyRequests.clear();
	}
}

void FindShortcuts(int n, struct weghtedGraph wGraph)
{
	std::vector<struct vertex> Q;
	for (int i = 0; i < n; i++)
	{
		struct vertex v;
		v.u = i;
		v.v = i;
		v.w = 0;
		Q.push_back(v);
	}
	
	
	while (!Q.empty())
	{
		std::multiset<struct vertex> q;
		for (int i = 0; i < Q.size(); i++)
		{
			struct vertex currentVertex = Q[i];
			int size = wGraph.sizes[currentVertex.v];
			std::vector<struct edge> neighbors = wGraph.graph[currentVertex.v];
			for (int j = 0; j < size; j++)
			{
				struct edge neighbor = neighbors[j];
				double neighbourWeight = neighbor.weight;
				//if (neighbourWeight > delta)
				//	break;
				struct vertex toInsert;
				toInsert.u = currentVertex.u;
				toInsert.v = neighbor.vertex;
				toInsert.w = currentVertex.w + neighbourWeight;
				q.insert(toInsert);
			}
		}
		
		std::multiset<struct vertex>::iterator it;
		struct vertex minVertex;
		minVertex.u = -1;
		minVertex.v = -1;
		minVertex.w = INT_MAX;
		
		Q.clear();
		for (it = q.begin(); it != q.end(); it++)
		{
			struct vertex currentVertex = *it;
			if (currentVertex.u != -1 && currentVertex.w <= delta && currentVertex.w < GetFromHash(currentVertex.u, currentVertex.v))
			{
				Q.push_back(currentVertex);
			}
			continue;
			
			
			//This should be a multiset distinction but..
			//if (currentVertex.u != minVertex.u || currentVertex.v != minVertex.v)
			//{
			//	if (minVertex.u != -1 && minVertex.w <= delta && minVertex.w < GetFromHash(minVertex.u, minVertex.v))
			//	{
			//		Q.push_back(minVertex);
			//	}
			//	
			//	minVertex = currentVertex;
			//}
			//else
			//{
			//	printf("Updating weight\n");
			//	if (minVertex.w > currentVertex.w)
			//		minVertex.w = currentVertex.w;
			//}
			//}
		}
		//if (minVertex.u != -1 && minVertex.w <= delta && minVertex.w < GetFromHash(minVertex.u, minVertex.v))
		//{
		//	Q.push_back(minVertex);
		//}
		
		for (int i = 0; i < Q.size(); i++)
		{
			struct vertex currentVertex = Q[i];
			PutToHash(currentVertex.u, currentVertex.v, currentVertex.w);
		}
	}
}

void DeltaSteppingImproved(int vertexCount, int startVertex, struct weghtedGraph wGraph, double* tent)
{
	InitializeHash(10);
	int totalHashSize = hashSize * hashSize;
	
	FindShortcuts(vertexCount, wGraph);
	
	printf("----\n");
	for (int i = 0; i < totalHashSize; i++)
	{
		std::vector<struct vertex> bucket = hash[i];
	
		for (int j = 0; j < bucket.size(); j++)
		{
			struct vertex currentVertex = bucket[j];
			struct edge newEdge;
			newEdge.vertex = currentVertex.v;
			newEdge.weight = currentVertex.w;
			//printf("Add edge %d -> %d and weight %f\n", currentVertex.u, newEdge.vertex, newEdge.weight);
			wGraph.graph[currentVertex.u].push_back(newEdge);
			wGraph.sizes[currentVertex.u]++;
		}
	}

	
	for (int i = 0; i < vertexCount; i++)
		tent[i] = INT_MAX;
	
	int bucketsCount = (double)1 / delta + 1;
	std::vector<std::set<int>> buckets (bucketsCount);
	
	
	struct request initialRequest;
	initialRequest.vertex = 0;
	initialRequest.weight = 0;
	relax(buckets, tent, initialRequest);
	
	int cyclingIndex = -1;
	while (true)
	{
		//printf("In while true\n");
		int minIndex = -1;
		for (int i = 0; i < bucketsCount; i++)
		{
			cyclingIndex++;
			int indexForCheck = cyclingIndex % bucketsCount;
			if (buckets[indexForCheck].size() > 0)
			{
				minIndex = indexForCheck;
				break;
			}
		}
		if (minIndex == -1)
		{
			break;
		}		
		
		//printf("Current cycling delta index = %d\n", cyclingIndex);
		
		std::set<int> currentSet = buckets[minIndex];
		//printf("Current bucket size = %d\n", currentSet.size());
		std::set<int>::iterator it;
		for (it = currentSet.begin(); it != currentSet.end(); it++)
		{
			int currentVertex = *it;
			std::vector<struct edge> neighbors = wGraph.graph[currentVertex];
			for (int i = 0; i < neighbors.size(); i++)
			{
				struct edge neighbor = neighbors[i];
				double desiredWeight = tent[currentVertex] + neighbor.weight;
				if (desiredWeight <= delta*(cyclingIndex + 1))
				{
					struct request req;
					req.vertex = neighbor.vertex;
					req.weight = desiredWeight;
					//printf("Relaxing request %d %f\n", req.vertex, req.weight);
					relax(buckets, tent, req);
				}
			}
		}
		for (it = currentSet.begin(); it != currentSet.end(); it++)
		{
			int currentVertex = *it;
			std::vector<struct edge> neighbors = wGraph.graph[currentVertex];
			for (int i = 0; i < neighbors.size(); i++)
			{
				struct edge neighbor = neighbors[i];
				double desiredWeight = tent[currentVertex] + neighbor.weight;
				if (desiredWeight > delta*(cyclingIndex + 1))
				{
					struct request req;
					req.vertex = neighbor.vertex;
					req.weight = desiredWeight;
					relax(buckets, tent, req);
				}
			}
		}
		
		buckets[minIndex].clear();
	}
}

int main(int argc, char *argv[])
{
    double time;
    int NumThreads = omp_get_num_procs();
	omp_set_nested(1);
        
    int startVertex = 0;
    int n = 20 * 1000;
	int maxDegree;
    
    if (sscanf (argv[1], "%d", &n) != 1) 
    {
        printf("error: use 'delta.exe n d' where n - vertices count and d - maximum degree\n");
    }
    if (sscanf (argv[2], "%d", &maxDegree) != 1) 
    {
        printf("error: use 'delta.exe n d' where n - vertices count and d - maximum degree\n");
    }
	
	std::cout.precision(10);
	
	delta = 0.1;
	
	
	
    time = omp_get_wtime();
	
    struct weghtedGraph wGraph = GenerateGraph(n, maxDegree);
    
    time = omp_get_wtime() - time;
    printf("Graph generating time: %f\n", time);
	
	
	
	
	//parallel code execution
    int* parentParallel = new int[n];
    double* pathWeightParallel = new double[n];      
       
    time = omp_get_wtime();
    
    FindShortestPathesParallel(n, startVertex, wGraph, parentParallel, pathWeightParallel);
    
    time = omp_get_wtime() - time;
    printf("%d threads Dijkstra calculating time: %f\n", NumThreads, time);
   
    
    
    
    
    
    //single thread code execution
    int* parentSingleThread = new int[n];
    double* pathWeightSingleThread = new double[n]; 
    
    time = omp_get_wtime();
    
    FindShortestPathes(n, startVertex, wGraph, parentSingleThread, pathWeightSingleThread);
    
    time = omp_get_wtime() - time;
    printf("Single thread Dijkstra calculating time: %f\n", time);
	
	
	
	
	double* pathWeightDeltaStepping = new double[n];
	
    time = omp_get_wtime();
	
	FindShortestPathesDeltaStepping(n, startVertex, wGraph, pathWeightDeltaStepping);
    
    time = omp_get_wtime() - time;
    printf("Delta-stepping calculation time: %f\n", time);
	
	
	
	
	
	double* pathWeightDeltaSteppingParallel = new double[n];
	
    time = omp_get_wtime();
	
	DeltaSteppingImproved(n, startVertex, wGraph, pathWeightDeltaSteppingParallel);
    
    time = omp_get_wtime() - time;
    printf("Delta-stepping improved parallel (%d threads) calculation time: %f\n", NumThreads, time);
	
	
	
	for (int i = 0; i < n; i++)
	{
		//std::cout << pathWeightParallel[i] << ' ' << pathWeightDeltaStepping[i] << ' ' << pathWeightDeltaSteppingParallel[i] << std::endl;
		if (pathWeightParallel[i] != pathWeightSingleThread[i])
		{
			printf("Smth wrong\n");
			return 1;
		}
		if (fabs(pathWeightParallel[i] - pathWeightDeltaStepping[i]) > 1e-5)
		{
			printf("pathWeightParallel[%d] = %f, but pathWeightDeltaStepping[%d] = %f\n", i, pathWeightParallel[i], i, pathWeightDeltaStepping[i]);
			return 1;
		}
		if (fabs(pathWeightParallel[i] - pathWeightDeltaSteppingParallel[i]) > 1e-5)
		{
			printf("pathWeightParallel[%d] = %f, but pathWeightDeltaSteppingParallel[%d] = %f\n", i, pathWeightParallel[i], i, pathWeightDeltaSteppingParallel[i]);
			return 1;
		}
	}
	
	printf("Results are equals!");
		
}