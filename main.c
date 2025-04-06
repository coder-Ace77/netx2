#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include <time.h>

#define MAX_NODES 5000
#define MAX_LINE_LENGTH 10000

typedef struct {
    int first;
    int second;
} Pair;

typedef struct {
    int** adj;          // Adjacency list
    int* degrees;       // Degree of each node
    int node_count;      // Total number of nodes
    Pair* edges;        // Array of edges
    int edge_count;     // Total number of edges
} Graph;

void initializeGraph(Graph* g) {
    g->adj = (int**)malloc(MAX_NODES * sizeof(int*));
    g->degrees = (int*)calloc(MAX_NODES, sizeof(int));
    for (int i = 0; i < MAX_NODES; i++) {
        g->adj[i] = (int*)malloc(MAX_NODES * sizeof(int));
    }
    g->edges = (Pair*)malloc(MAX_NODES * MAX_NODES * sizeof(Pair));
    g->edge_count = 0;
    g->node_count = 0;
}

void freeGraph(Graph* g) {
    for (int i = 0; i < MAX_NODES; i++) {
        free(g->adj[i]);
    }
    free(g->adj);
    free(g->degrees);
    free(g->edges);
}

void readGraph(Graph* g, const char* filename) {
    printf("Reading graph from %s\n", filename);
    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    char line[MAX_LINE_LENGTH];
    int max_node = -1;

    printf("Reading graph\n");
    // Read the graph from the file

    while (fgets(line, sizeof(line), file)){
        printf("Reading line: %s", line);
        if (line[0] == '#' || line[0] == '\n') continue;

        int u;
        char* ptr = line;
        u = strtol(ptr, &ptr, 10);

        if (u > max_node) max_node = u;

        int v;
        while (*ptr != '\n'){
            v = strtol(ptr, &ptr, 10);
            if (u == v) continue;
            // Add to adjacency matrix
            g->adj[u][g->degrees[u]++] = v;
            g->adj[v][g->degrees[v]++] = u;

            // Add to edges list (only once per pair)
            if (u < v) {
                g->edges[g->edge_count].first = u;
                g->edges[g->edge_count].second = v;
                g->edge_count++;
            }

            if (v > max_node) max_node = v;
        }
    }

    g->node_count = max_node + 1;
    fclose(file);
}

void addEdge(Graph* g, int u, int v) {
    g->adj[u][g->degrees[u]++] = v;
    g->adj[v][g->degrees[v]++] = u;
}

void removeEdge(Graph* g, int u, int v) {
    // Remove v from u's adjacency list
    for (int i = 0; i < g->degrees[u]; i++) {
        if (g->adj[u][i] == v) {
            g->adj[u][i] = g->adj[u][g->degrees[u] - 1];
            g->degrees[u]--;
            break;
        }
    }

    // Remove u from v's adjacency list
    for (int i = 0; i < g->degrees[v]; i++) {
        if (g->adj[v][i] == u) {
            g->adj[v][i] = g->adj[v][g->degrees[v] - 1];
            g->degrees[v]--;
            break;
        }
    }
}

Pair computeCentralities(Graph* g){
    double* betweenness = (double*)calloc(g->node_count, sizeof(double));
    long long* stress = (long long*)calloc(g->node_count, sizeof(long long));

    for (int s = 0; s < g->node_count; s++) {
        if (g->degrees[s] == 0) continue;

        int* pred[MAX_NODES];
        int pred_count[MAX_NODES] = {0};
        for (int i = 0; i < g->node_count; i++) {
            pred[i] = (int*)malloc(g->degrees[i] * sizeof(int));
        }

        int* dist = (int*)malloc(g->node_count * sizeof(int));
        long long* sigma = (long long*)calloc(g->node_count, sizeof(long long));
        int* stack = (int*)malloc(g->node_count * sizeof(int));
        int stack_top = -1;
        int* queue = (int*)malloc(g->node_count * sizeof(int));
        int queue_front = 0, queue_rear = 0;

        // Initialize distances and sigma
        for (int i = 0; i < g->node_count; i++) {
            dist[i] = -1;
        }
        dist[s] = 0;
        sigma[s] = 1;
        queue[queue_rear++] = s;

        // BFS phase
        while (queue_front < queue_rear) {
            int v = queue[queue_front++];
            stack[++stack_top] = v;

            for (int i = 0; i < g->degrees[v]; i++) {
                int w = g->adj[v][i];
                if (dist[w] < 0) {
                    dist[w] = dist[v] + 1;
                    queue[queue_rear++] = w;
                }
                if (dist[w] == dist[v] + 1) {
                    sigma[w] += sigma[v];
                    pred[w][pred_count[w]++] = v;
                }
            }
        }

        // Accumulation phase
        double* delta = (double*)calloc(g->node_count, sizeof(double));
        while (stack_top >= 0) {
            int w = stack[stack_top--];
            for (int i = 0; i < pred_count[w]; i++) {
                int v = pred[w][i];
                double contrib = ((double)sigma[v] * (1.0 + delta[w])) / sigma[w];
                delta[v] += contrib;

                if (v != s) {
                    stress[v] += sigma[v];
                }
            }
            if (w != s) {
                betweenness[w] += delta[w] / 2.0;
            }
        }

        for (int i = 0; i < g->node_count; i++) {
            free(pred[i]);
        }
        free(dist);
        free(sigma);
        free(stack);
        free(queue);
        free(delta);
    }

    // Find maximum values
    double max_betweenness = 0.0;
    long long max_stress = 0;
    for (int i = 0; i < g->node_count; i++) {
        if (betweenness[i] > max_betweenness) {
            max_betweenness = betweenness[i];
        }
        if (stress[i] > max_stress) {
            max_stress = stress[i];
        }
    }

    free(betweenness);
    free(stress);

    Pair result = {max_betweenness, max_stress};
    return result;
}

void findMissingEdges(Graph* g, Pair** missing_edges, int* missing_count) {
    bool** existing = (bool**)calloc(g->node_count, sizeof(bool*));
    for (int i = 0; i < g->node_count; i++) {
        existing[i] = (bool*)calloc(g->node_count, sizeof(bool));
    }

    // Mark existing edges
    for (int i = 0; i < g->edge_count; i++) {
        int u = g->edges[i].first;
        int v = g->edges[i].second;
        existing[u][v] = true;
        existing[v][u] = true;
    }

    // Count missing edges
    *missing_count = 0;
    for (int u = 0; u < g->node_count; u++) {
        for (int v = u + 1; v < g->node_count; v++) {
            if (!existing[u][v]) {
                (*missing_count)++;
            }
        }
    }

    // Allocate and store missing edges
    *missing_edges = (Pair*)malloc(*missing_count * sizeof(Pair));
    int index = 0;
    for (int u = 0; u < g->node_count; u++) {
        for (int v = u + 1; v < g->node_count; v++) {
            if (!existing[u][v]) {
                (*missing_edges)[index].first = u;
                (*missing_edges)[index].second = v;
                index++;
            }
        }
    }

    for (int i = 0; i < g->node_count; i++) {
        free(existing[i]);
    }
    free(existing);
}

void solve(const char* filename){
    Graph g;
    initializeGraph(&g);
    readGraph(&g, filename);

    clock_t start = clock();
    Pair initial = computeCentralities(&g);
    printf("Initial betweenness: %.2f, stress: %lld\n", initial.first, initial.second);

    Pair* missing_edges;
    int missing_count;
    findMissingEdges(&g, &missing_edges, &missing_count);
    printf("Total missing edges: %d\n", missing_count);

    double min_betweenness = initial.first;
    long long min_stress = initial.second;
    Pair best_betweenness_edge = {-1, -1};
    Pair best_stress_edge = {-1, -1};

    for (int i = 0; i < missing_count; i++) {
        Pair edge = missing_edges[i];
        addEdge(&g, edge.first, edge.second);

        Pair current = computeCentralities(&g);

        if (current.first < min_betweenness) {
            min_betweenness = current.first;
            best_betweenness_edge = edge;
        }

        if (current.second < min_stress) {
            min_stress = current.second;
            best_stress_edge = edge;
        }

        removeEdge(&g, edge.first, edge.second);

        if (i % 1000 == 0) {
            printf("Processed %d edges...\n", i);
        }
    }

    printf("\nFinal results:\n");
    printf("Minimum betweenness: %.2f at edge (%d, %d)\n", 
           min_betweenness, best_betweenness_edge.first, best_betweenness_edge.second);
    printf("Minimum stress: %lld at edge (%d, %d)\n", 
           min_stress, best_stress_edge.first, best_stress_edge.second);

    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Total execution time: %.2f seconds\n", time_spent);

    free(missing_edges);
    freeGraph(&g);
}

void main() {
    printf("Graph Processing Program\n");
    solve("graph1.adjlist");
}