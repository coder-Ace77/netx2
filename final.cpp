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
#include <unordered_map>
#include <unordered_set>

using namespace std;
using namespace std::chrono;

class Graph {
public:
    vector<vector<int>> adj;
    vector<pair<int, int>> edges;  
    void readGraph(const string& filename){
        ifstream file(filename);
        if (!file.is_open()) {
            cerr << "Error opening file: " << filename << endl;
            return;
        }

        edges.clear();
        string line;
        int maxNode = 0;
        set<pair<int,int>> edgeSet;
        while (getline(file, line)){
            if (line.empty() || line[0] == '#') continue;
    
            istringstream iss(line);
            int u;
            if (!(iss >> u)) continue;
    
            int v;
            while (iss >> v) {
                if (u == v) continue;
                int b = max(u, v);
                maxNode = max(maxNode, b);
                if(u>v) swap(u, v);
                edgeSet.insert({u,v});
            }
        }
    
        adj.clear();
        adj.resize(maxNode + 1);
        for (auto [u, v] : edgeSet){
            edges.push_back({u, v});
            adj[u].push_back(v);
            adj[v].push_back(u);
        }
        file.close();
    }
    
    void addEdge(int u, int v){
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    
    vector<pair<int, int>> findMissingEdges(){
        vector<pair<int, int>> missingEdges;
        set<pair<int,int>> existingEdges(edges.begin(), edges.end());
        for (int u = 0; u < adj.size(); ++u) {
            for (int v = u+1; v < adj.size(); ++v){
                if(existingEdges.count({u,v}) == 0){
                    missingEdges.push_back({u,v});
                }
            }
        }
        return missingEdges;
    }
    
    vector<vector<int>> pred;

    inline pair<double, long long> compute(bool test=false){
        int V = adj.size();
        if(pred.size()<V){
            pred.resize(V,vector<int> (V));
        }
        double betweenness[V]={0.0};
        long long stress[V]={0};
        int dist[V]={-1};
        long long sigma[V]={0}; 
        int predIND[V]={0};
        stack<int> S;
        queue<int> Q;
        double max_bet = 0.0;
        long long max_stress = 0;
        for (int s=0;s<V;++s){
            fill(predIND,predIND+V,0);
            fill(dist,dist+V,-1);
            fill(sigma,sigma+V,0ll); 
            dist[s]=0;
            sigma[s]=1;
            Q.push(s);
            while(!Q.empty()){
                int v = Q.front();
                Q.pop();
                S.push(v);
                for(int i=0;i<adj[v].size();i++){
                    int w = adj[v][i];
                    if (dist[w]<0){
                        dist[w]=dist[v]+1;
                        Q.push(w);
                    }
                    if (dist[w]==dist[v]+1){
                        sigma[w]+=sigma[v];
                        pred[w][predIND[w]]=v;
                        ++predIND[w];
                    }
                }
            }    
            double delta[V]={0.0};
            long long delta2[V]={0};
            while (!S.empty()){
                int w = S.top();
                S.pop();
                for (int i=0;i<predIND[w];i++){
                    int v = pred[w][i];
                    double contrib = ((double)sigma[v]*(1.0 + delta[w]))/sigma[w];
                    delta[v] += contrib;
                    delta2[v] += (delta2[w]+1);
                }
                if (w != s) {
                    betweenness[w]+=delta[w];
                    stress[w]+=delta2[w]*sigma[w];
                }              
            }
        }
        for (int i=0;i<V;++i){
            betweenness[i]/=2.0;
            stress[i]/=2;
            if (betweenness[i] > max_bet){
                max_bet = betweenness[i];
            }
            if(stress[i]>max_stress){
                max_stress = stress[i];
            }
        }
        if(test){
            for (int i=0;i<V;++i)cout << stress[i] << " ";
        }
        return {max_bet, max_stress};
    }
};

void solve_parallel(const string& filename,bool metrics = false , bool test = false){
    Graph g;
    g.readGraph(filename);
    auto ans = make_pair(1e9, (int)1e9);
    ans = g.compute();
    auto missingEdges = g.findMissingEdges();
    pair<int, int> ed1 = {-1, -1};
    pair<int, int> ed2 = {-1, -1};
    double min_betweenness = ans.first;
    long long min_stress = ans.second;

    #pragma omp parallel
    {
        Graph local_g = g;
        double local_min_betweenness = ans.first;
        long long local_min_stress = ans.second;
        pair<int, int> local_ed1 = {-1, -1};
        pair<int, int> local_ed2 = {-1, -1};

        #pragma omp for schedule(dynamic)
        for (int i = 0; i < (int)missingEdges.size(); ++i){
            auto x = missingEdges[i];
            local_g.addEdge(x.first, x.second);
            auto temp = local_g.compute();

            if (temp.first < local_min_betweenness) {
                local_min_betweenness = temp.first;
                local_ed1 = x;
            }
            if (temp.second < local_min_stress) {
                local_min_stress = temp.second;
                local_ed2 = x;
            }

            local_g.adj[x.first].pop_back();
            local_g.adj[x.second].pop_back();
            #pragma omp critical
            {
                if (i % 1000 == 0) {
                    cout<<"#"<<flush;
                }
            }
        }

        #pragma omp critical
        {
            if (local_min_betweenness < min_betweenness) {
                min_betweenness = local_min_betweenness;
                ed1 = local_ed1;
            }
            if (local_min_stress < min_stress) {
                min_stress = local_min_stress;
                ed2 = local_ed2;
            }
        }
    }
    cout<<endl;
    cout << "\nFinal Results:" << endl;
    cout << "Minimum Betweenness Centrality: " << min_betweenness << endl;
    cout << "Minimum Stress Centrality: " << min_stress << endl;
    cout << "BETWEENNESS minimized at edge: (" << ed1.first << ", " << ed1.second << ")" << endl;
    cout << "STRESS minimized at edge: (" << ed2.first << ", " << ed2.second << ")" << endl;

}

void solve_sequential(const string& filename){
    Graph g;
    g.readGraph(filename);
    auto ans = make_pair(1e9,(int)1e9);
    ans = g.compute();
    cout << "Initial betweenness centrality, stress centrality : " << ans.first<<" "<<ans.second << endl;
    auto missingEdges = g.findMissingEdges();
    cout << "Total missing edges to process: " << missingEdges.size() << endl;

    pair<int, int> ed1 = {-1, -1};
    pair<int, int> ed2 = {-1, -1};
    double min_betweenness = ans.first;
    long long min_stress = ans.second;
    int cnt = 0;


    for (size_t i = 0; i < missingEdges.size(); ++i) {
        auto x = missingEdges[i];
        g.addEdge(x.first, x.second);
        auto temp = g.compute();

        if(temp.first < min_betweenness){
            min_betweenness = temp.first;
            ed1 = x;
        }
        if (temp.second < min_stress){
            min_stress = temp.second;
            ed2 = x;
        }
        g.adj[x.first].pop_back();
        g.adj[x.second].pop_back();
        ans = {min_betweenness, min_stress};
        ++cnt;
        if(cnt%100==0)cout<<cnt<<" done\n";
    }

    cout << "\nFinal Results:" << endl;

    cout << "Minimum Betweenness Centrality: " << ans.first << endl;
    cout << "Minimum Stress Centrality: " << ans.second << endl;

    cout << "BETWEENNESS minimized at edge: (" << ed1.first << ", " << ed1.second << ")" << endl;
    cout << "STRESS minimized at edge: (" << ed2.first << ", " << ed2.second << ")" << endl;
}

int main(int argc, char* argv[]){
    string filename;
    if (argc < 2){
        cout<<"ENTER THE FILENAME: ";
        cin>>filename;
    }else{
        filename = argv[1];
    }
    solve_parallel(filename);
    return 0;   
}