#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <fstream>
#include <sstream>
#include <set>
#include <utility>
#include <chrono>
#include <omp.h>

using namespace std;
using namespace std::chrono;

class Graph {
public:
    vector<vector<int>> adj;
    vector<pair<int, int>> edges;  
    
    void readGraph(const string& filename) {
        ifstream file(filename);
        if (!file.is_open()) {
            cerr << "Error opening file: " << filename << endl;
            return;
        }
        string line;
        while (getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;
            break;
        }
        edges.clear();
        while(getline(file, line)) {
            if (line.empty()) continue;
            istringstream iss(line);
            int v;
            adj.push_back({});
            int u;
            iss >> u;
            while (iss >> v) {
                adj.back().push_back(v);
                if (u < v) {
                    edges.emplace_back(u, v);
                }
            }
        }
        file.close();
    }

    void addEdge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    
    vector<pair<int, int>> findMissingEdges() {
        vector<pair<int, int>> missingEdges;
        set<pair<int,int>> existingEdges(edges.begin(), edges.end());
        cout << "Existing edges: " << existingEdges.size() << endl;
        for (int u = 0; u < adj.size(); ++u) {
            for (int v = u+1; v < adj.size(); ++v) {
                if(existingEdges.count({u,v}) == 0) {
                    missingEdges.push_back({u,v});
                }
            }
        }
        return missingEdges;
    }

    pair<double, long long> compute() {
        int V = adj.size();
        vector<double> betweenness(V, 0.0);
        vector<long long> stress(V, 0);
    
        for (int s = 0; s < V; ++s) {
            vector<vector<int>> pred(V);    // predecessors on shortest paths
            vector<int> dist(V, -1);        // distance from source s
            vector<long long> sigma(V, 0);  // number of shortest paths from s
            stack<int> S;                  // stack of vertices in order of non-decreasing distance
            queue<int> Q;                   // for BFS
    
            dist[s] = 0;
            sigma[s] = 1;
            Q.push(s);
    
            // BFS phase: compute shortest paths and predecessors
            while (!Q.empty()) {
                int v = Q.front();
                Q.pop();
                S.push(v);
    
                for (int w : adj[v]) {
                    if (dist[w] < 0) {
                        dist[w] = dist[v] + 1;
                        Q.push(w);
                    }
                    if (dist[w] == dist[v] + 1) {
                        sigma[w] += sigma[v];
                        pred[w].push_back(v);
                    }
                }
            }
    
            // Accumulation phase: compute dependencies and contribute to stress
            vector<double> delta(V, 0.0);
            while (!S.empty()) {
                int w = S.top();
                S.pop();
                for (int v : pred[w]) {
                    double contrib = ((double)sigma[v] * (1.0 + delta[w])) / sigma[w];
                    delta[v] += contrib;
    
                    // Update stress centrality for the predecessor v,
                    // since v is used to forward shortest paths to w.
                    if (v != s) {
                        stress[v] += sigma[v];
                    }
                }
                if (w != s) {
                    betweenness[w] += delta[w] / 2.0;  // adjustment for undirected graphs
                }
            }
        }
    
        double max_bet = 0.0;
        long long max_stress = 0;
        for (int i = 0; i < V; ++i) {
            max_bet = max(max_bet, betweenness[i]);
            max_stress = max(max_stress, stress[i]);
        }
    
        return {max_bet, max_stress};
    }
};

void solve(string filename){
    Graph g;
    g.readGraph("graph1.adjlist");
    
    auto totalStart = high_resolution_clock::now();
    auto ans = g.compute();

    int num_procs = omp_get_num_procs();
    cout << "Number of processors: " << num_procs << endl;
    
    auto missingEdges = g.findMissingEdges();
    cout << "Total missing edges to process: " << missingEdges.size() << endl;
    
    pair<int,int> ed1 = {-1,-1};
    pair<int,int> ed2 = {-1,-1};
    
    double min_betweenness = ans.first;
    long long min_stress = ans.second;
    
    #pragma omp parallel
    {
        Graph local_g = g;        
        pair<int,int> local_ed1 = {-1,-1};
        pair<int,int> local_ed2 = {-1,-1};
        double local_min_bet = min_betweenness;
        long long local_min_stress = min_stress;
        
        #pragma omp for schedule(dynamic) nowait
        for (size_t i = 0; i < missingEdges.size(); ++i) {
            auto x = missingEdges[i];
            
            // Add the edge to the local graph copy
            local_g.addEdge(x.first, x.second);
            
            // Compute metrics with the new edge
            auto temp = local_g.compute();
            
            // Check for new minimum betweenness
            if (temp.first < local_min_bet) {
                local_min_bet = temp.first;
                local_ed1 = x;
            }
            
            // Check for new minimum stress
            if (temp.second < local_min_stress) {
                local_min_stress = temp.second;
                local_ed2 = x;
            }
            
            // Remove the edge from the local graph copy
            local_g.adj[x.first].pop_back();
            local_g.adj[x.second].pop_back();
            
            if(i % 1000 == 0){
                #pragma omp critical
                {
                    cout << i << " done by thread " << omp_get_thread_num() << "\n" << flush;
                }
            }
        }
        
        // Update global minimums with thread-local results
        #pragma omp critical
        {
            if (local_min_bet < min_betweenness) {
                min_betweenness = local_min_bet;
                ed1 = local_ed1;
            }
            if (local_min_stress < min_stress) {
                min_stress = local_min_stress;
                ed2 = local_ed2;
            }
        }
    }
    
    ans = {min_betweenness, min_stress};
    
    cout << ans.first << " " << ans.second << endl;
    cout << "BETWEENESS minimised at " << ed1.first << " " << ed1.second << endl;
    cout << "STRESS minimised at " << ed2.first << " " << ed2.second << endl;
    
    auto totalEnd = high_resolution_clock::now();
    auto totalDuration = duration_cast<seconds>(totalEnd - totalStart);
    cout << "Total execution time: " << totalDuration.count() << " s" << endl;
}

int main(){
    for(int i=1;i<=5;i++){
        string ss = "graph"+to_string(i)+".adjlist";
        solve(ss);
    }
    return 0;   
}