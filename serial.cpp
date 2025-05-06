#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <climits>

using namespace std;

struct BipartiteGraph {
    int num_u;
    int num_v;
    int num_edges;
    vector<int> u_offset;
    vector<int> u_edges;
};

struct VertexDegree {
    int u;
    int degree;
};

struct VVertexInfo {
    int v;
    int degree;
    int rank;
};

struct EdgeButterfly {
    int u;
    int v;
    long long butterflies;
};

BipartiteGraph load_graph(const string& filename) {
    BipartiteGraph graph;
    ifstream file(filename);
    string line;
    if (!file.is_open()) {
        cerr << "Error: Could not open " << filename << endl;
        exit(1);
    }
    getline(file, line);
    if (line != "AdjacencyGraph") {
        cerr << "Error: Invalid header, expected 'AdjacencyGraph', got '" << line << "'\n";
        exit(1);
    }
    getline(file, line);
    graph.num_u = stoi(line);
    getline(file, line);
    graph.num_edges = stoi(line);
    graph.u_offset.resize(graph.num_u + 1, 0);
    graph.u_edges.clear();
    vector<int> all_v_indices;
    int actual_edges = 0;
    for (int u = 0; u < graph.num_u; ++u) {
        if (!getline(file, line)) {
            cerr << "Error: Incomplete graph data at vertex " << u << endl;
            exit(1);
        }
        stringstream ss(line);
        int v, edge_count = 0;
        while (ss >> v) {
            graph.u_edges.push_back(v);
            all_v_indices.push_back(v);
            edge_count++;
        }
        actual_edges += edge_count;
        graph.u_offset[u + 1] = graph.u_edges.size();
    }
    if (actual_edges != graph.num_edges) {
        cerr << "Error: Expected " << graph.num_edges << " edges, but loaded " << actual_edges << endl;
        exit(1);
    }
    if (all_v_indices.empty()) {
        graph.num_v = 0;
    } else {
        graph.num_v = *max_element(all_v_indices.begin(), all_v_indices.end()) + 1;
    }
    file.close();
    return graph;
}

void preprocess_graph(BipartiteGraph& graph) {
    vector<int> v_degrees(graph.num_v, 0);
    for (int e : graph.u_edges) {
        v_degrees[e]++;
    }
    vector<VVertexInfo> v_info(graph.num_v);
    for (int v = 0; v < graph.num_v; ++v) {
        v_info[v] = {v, v_degrees[v], 0};
    }
    sort(v_info.begin(), v_info.end(), 
         [](const VVertexInfo& a, const VVertexInfo& b) {
             return a.degree > b.degree;
         });
    for (int i = 0; i < graph.num_v; ++i) {
        v_info[i].rank = i;
    }
    vector<int> v_rank(graph.num_v);
    for (const auto& vi : v_info) {
        v_rank[vi.v] = vi.rank;
    }
    vector<int> new_u_edges;
    vector<int> new_u_offset(graph.num_u + 1, 0);
    new_u_edges.reserve(graph.u_edges.size());
    for (int u = 0; u < graph.num_u; ++u) {
        int start = graph.u_offset[u];
        int end = graph.u_offset[u + 1];
        vector<pair<int, int>> ranked_neighbors;
        for (int i = start; i < end; ++i) {
            int v = graph.u_edges[i];
            ranked_neighbors.emplace_back(v_rank[v], v);
        }
        sort(ranked_neighbors.begin(), ranked_neighbors.end());
        for (const auto& rn : ranked_neighbors) {
            new_u_edges.push_back(rn.second);
        }
        new_u_offset[u + 1] = new_u_edges.size();
    }
    graph.u_edges = new_u_edges;
    graph.u_offset = new_u_offset;
}

vector<EdgeButterfly> count_butterflies(const BipartiteGraph& graph) {
    vector<EdgeButterfly> butterflies;
    vector<vector<int>> v_to_u(graph.num_v);
    for (int u = 0; u < graph.num_u; ++u) {
        int start = graph.u_offset[u];
        int end = graph.u_offset[u + 1];
        for (int i = start; i < end; ++i) {
            v_to_u[graph.u_edges[i]].push_back(u);
        }
    }
    for (int u = 0; u < graph.num_u; ++u) {
        int start = graph.u_offset[u];
        int end = graph.u_offset[u + 1];
        for (int i = start; i < end; ++i) {
            int v = graph.u_edges[i];
            long long butterfly_count = 0;
            for (int j = i + 1; j < end; ++j) {
                int vp = graph.u_edges[j];
                for (int up : v_to_u[vp]) {
                    if (up != u) {
                        butterfly_count++;
                    }
                }
            }
            butterflies.push_back({u, v, butterfly_count});
        }
    }
    return butterflies;
}

BipartiteGraph peel_edges(BipartiteGraph& graph, const vector<EdgeButterfly>& butterflies, long long min_butterflies) {
    vector<pair<int, int>> edges_to_remove;
    for (const auto& eb : butterflies) {
        if (eb.butterflies == min_butterflies) {
            edges_to_remove.emplace_back(eb.u, eb.v);
        }
    }
    vector<vector<int>> new_adj_list(graph.num_u);
    for (int u = 0; u < graph.num_u; ++u) {
        int start = graph.u_offset[u];
        int end = graph.u_offset[u + 1];
        for (int i = start; i < end; ++i) {
            int v = graph.u_edges[i];
            bool keep = true;
            for (const auto& er : edges_to_remove) {
                if (er.first == u && er.second == v) {
                    keep = false;
                    break;
                }
            }
            if (keep) {
                new_adj_list[u].push_back(v);
            }
        }
    }
    BipartiteGraph new_graph;
    new_graph.num_u = graph.num_u;
    new_graph.num_v = graph.num_v;
    new_graph.u_offset.resize(graph.num_u + 1, 0);
    new_graph.u_edges.clear();
    int new_edges = 0;
    for (int u = 0; u < graph.num_u; ++u) {
        new_graph.u_edges.insert(new_graph.u_edges.end(), new_adj_list[u].begin(), new_adj_list[u].end());
        new_edges += new_adj_list[u].size();
        new_graph.u_offset[u + 1] = new_graph.u_edges.size();
    }
    new_graph.num_edges = new_edges;
    return new_graph;
}

int main() {
    auto start_time = chrono::high_resolution_clock::now();
    cout << "====================================\n";
    cout << "Serial Butterfly Counting with Peeling\n";
    cout << "====================================\n";
    
    string filename = "/home/storage/graph.txt";
    BipartiteGraph graph = load_graph(filename);
    
    int max_iterations = 10;
    for (int iter = 0; iter < max_iterations; ++iter) {
        if (graph.num_edges == 0) {
            cout << "\nIteration " << iter << ": No edges remaining\n";
            break;
        }
        preprocess_graph(graph);
        
        cout << "\nIteration " << iter << " Statistics:\n";
        cout << "  U Vertices: " << graph.num_u << "\n";
        cout << "  V Vertices: " << graph.num_v << "\n";
        cout << "  Edges: " << graph.num_edges << "\n";
        
        vector<EdgeButterfly> butterflies = count_butterflies(graph);
        long long min_butterflies = LLONG_MAX;
        for (const auto& eb : butterflies) {
            if (eb.butterflies > 0 && eb.butterflies < min_butterflies) {
                min_butterflies = eb.butterflies;
            }
        }
        if (min_butterflies == LLONG_MAX) {
            cout << "\nIteration " << iter << ": No butterflies remaining\n";
            break;
        }
        
        sort(butterflies.begin(), butterflies.end(), 
             [](const EdgeButterfly& a, const EdgeButterfly& b) {
                 return a.u == b.u ? a.v < b.v : a.u < b.u;
             });
        long long total_butterflies = 0;
        for (const auto& eb : butterflies) {
            total_butterflies += eb.butterflies;
        }
        cout << "\nIteration " << iter << " Butterfly Counts:\n";
        cout << "----------------------------------------\n";
        cout << "| " << left << setw(6) << "U" << " | " << setw(6) << "V" << " | " << setw(12) << "Butterflies" << " |\n";
        cout << "----------------------------------------\n";
        for (const auto& eb : butterflies) {
            if (eb.butterflies > 0) {
                cout << "| " << left << setw(6) << eb.u << " | " << setw(6) << eb.v << " | " << setw(12) << eb.butterflies << " |\n";
            }
        }
        cout << "----------------------------------------\n";
        cout << "  Total Butterflies: " << total_butterflies << "\n";
        cout << "  Minimum Butterfly Count: " << min_butterflies << "\n";
        
        graph = peel_edges(graph, butterflies, min_butterflies);
    }
    
    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
    cout << "\nFinal Summary:\n";
    cout << "  Remaining Edges: " << graph.num_edges << "\n";
    cout << "  Execution Time: " << fixed << setprecision(3) << duration.count() / 1000.0 << " seconds\n";
    cout << "====================================\n";
    
    return 0;
}