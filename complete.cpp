#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <mpi.h>
#include <algorithm>
#include <iomanip>
#include <omp.h>
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

struct LocalGraph {
    vector<int> local_u;
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
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    getline(file, line);
    if (line != "AdjacencyGraph") {
        cerr << "Error: Invalid header, expected 'AdjacencyGraph', got '" << line << "'\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
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
            MPI_Abort(MPI_COMM_WORLD, 1);
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
        MPI_Abort(MPI_COMM_WORLD, 1);
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

LocalGraph partition_graph(const BipartiteGraph& graph, int rank, int size) {
    LocalGraph local;
    vector<VertexDegree> degrees(graph.num_u);
    for (int u = 0; u < graph.num_u; ++u) {
        int degree = graph.u_offset[u + 1] - graph.u_offset[u];
        degrees[u] = {u, degree};
    }
    sort(degrees.begin(), degrees.end(), 
         [](const VertexDegree& a, const VertexDegree& b) {
             return a.degree > b.degree;
         });
    vector<int> edge_counts(size, 0);
    vector<vector<int>> node_vertices(size);
    int target_edges = graph.num_edges / size;
    for (const auto& vd : degrees) {
        int min_node = 0;
        int min_edges = edge_counts[0];
        for (int i = 1; i < size; ++i) {
            if (edge_counts[i] < min_edges) {
                min_node = i;
                min_edges = edge_counts[i];
            }
        }
        node_vertices[min_node].push_back(vd.u);
        edge_counts[min_node] += vd.degree;
    }
    local.local_u = node_vertices[rank];
    sort(local.local_u.begin(), local.local_u.end());
    local.u_offset.push_back(0);
    int edge_count = 0;
    for (int u : local.local_u) {
        int start = graph.u_offset[u];
        int end = graph.u_offset[u + 1];
        for (int i = start; i < end; ++i) {
            local.u_edges.push_back(graph.u_edges[i]);
            edge_count++;
        }
        local.u_offset.push_back(local.u_edges.size());
    }
    return local;
}

vector<EdgeButterfly> count_butterflies(const BipartiteGraph& graph, const LocalGraph& local, int rank) {
    vector<EdgeButterfly> local_butterflies;
    vector<vector<int>> v_to_u(graph.num_v);
    for (int u_idx = 0; u_idx < local.local_u.size(); ++u_idx) {
        int u = local.local_u[u_idx];
        int start = local.u_offset[u_idx];
        int end = local.u_offset[u_idx + 1];
        for (int i = start; i < end; ++i) {
            v_to_u[local.u_edges[i]].push_back(u);
        }
    }
    #pragma omp parallel
    {
        vector<EdgeButterfly> thread_butterflies;
        #pragma omp for schedule(dynamic)
        for (int u_idx = 0; u_idx < local.local_u.size(); ++u_idx) {
            int u = local.local_u[u_idx];
            int start = local.u_offset[u_idx];
            int end = local.u_offset[u_idx + 1];
            for (int i = start; i < end; ++i) {
                int v = local.u_edges[i];
                long long butterfly_count = 0;
                for (int j = i + 1; j < end; ++j) {
                    int vp = local.u_edges[j];
                    for (int up : v_to_u[vp]) {
                        if (up != u) {
                            butterfly_count++;
                        }
                    }
                }
                thread_butterflies.push_back({u, v, butterfly_count});
            }
        }
        #pragma omp critical
        local_butterflies.insert(local_butterflies.end(), thread_butterflies.begin(), thread_butterflies.end());
    }
    return local_butterflies;
}

BipartiteGraph peel_edges(BipartiteGraph& graph, const vector<EdgeButterfly>& butterflies, long long min_butterflies, int rank, int size) {
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

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    double start_time = MPI_Wtime();
    cout << flush;
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        cout << "====================================\n";
        cout << "Butterfly Counting with Peeling\n";
        cout << "====================================\n";
        cout << "MPI Processes: " << size << "\n";
    }
    BipartiteGraph graph;
    if (rank == 0) {
        string filename = "/home/storage/graph.txt";
        graph = load_graph(filename);
    }
    MPI_Bcast(&graph.num_u, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&graph.num_v, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&graph.num_edges, 1, MPI_INT, 0, MPI_COMM_WORLD);
    vector<int> u_offset(graph.num_u + 1);
    vector<int> u_edges(graph.num_edges);
    if (rank == 0) {
        u_offset = graph.u_offset;
        u_edges = graph.u_edges;
    }
    MPI_Bcast(u_offset.data(), graph.num_u + 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(u_edges.data(), graph.num_edges, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        graph.u_offset = u_offset;
        graph.u_edges = u_edges;
    }
    int max_iterations = 10;
    for (int iter = 0; iter < max_iterations; ++iter) {
        if (graph.num_edges == 0) {
            if (rank == 0) {
                cout << "\nIteration " << iter << ": No edges remaining\n";
            }
            break;
        }
        preprocess_graph(graph);
        LocalGraph local_graph = partition_graph(graph, rank, size);
        int local_edges = local_graph.u_edges.size();
        int local_vertices = local_graph.local_u.size();
        vector<int> all_local_edges(size);
        vector<int> all_local_vertices(size);
        MPI_Gather(&local_edges, 1, MPI_INT, all_local_edges.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gather(&local_vertices, 1, MPI_INT, all_local_vertices.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            cout << "\nIteration " << iter << " Statistics:\n";
            cout << "  U Vertices: " << graph.num_u << "\n";
            cout << "  V Vertices: " << graph.num_v << "\n";
            cout << "  Edges: " << graph.num_edges << "\n";
            cout << "  Partitioning:\n";
            for (int i = 0; i < size; ++i) {
                cout << "    Rank " << i << ": " << all_local_vertices[i] << " U vertices, " << all_local_edges[i] << " edges\n";
            }
        }
        vector<EdgeButterfly> local_butterflies = count_butterflies(graph, local_graph, rank);
        vector<EdgeButterfly> all_butterflies;
        long long min_butterflies = LLONG_MAX;
        if (rank == 0) {
            all_butterflies = local_butterflies;
            for (int r = 1; r < size; ++r) {
                int count;
                MPI_Recv(&count, 1, MPI_INT, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                vector<int> recv_u(count);
                vector<int> recv_v(count);
                vector<long long> recv_bf(count);
                MPI_Recv(recv_u.data(), count, MPI_INT, r, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(recv_v.data(), count, MPI_INT, r, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(recv_bf.data(), count, MPI_LONG_LONG, r, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (int i = 0; i < count; ++i) {
                    all_butterflies.push_back({recv_u[i], recv_v[i], recv_bf[i]});
                }
            }
            for (const auto& eb : all_butterflies) {
                if (eb.butterflies > 0 && eb.butterflies < min_butterflies) {
                    min_butterflies = eb.butterflies;
                }
            }
        } else {
            int count = local_butterflies.size();
            MPI_Send(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            vector<int> send_u(count);
            vector<int> send_v(count);
            vector<long long> send_bf(count);
            for (int i = 0; i < count; ++i) {
                send_u[i] = local_butterflies[i].u;
                send_v[i] = local_butterflies[i].v;
                send_bf[i] = local_butterflies[i].butterflies;
            }
            MPI_Send(send_u.data(), count, MPI_INT, 0, 1, MPI_COMM_WORLD);
            MPI_Send(send_v.data(), count, MPI_INT, 0, 2, MPI_COMM_WORLD);
            MPI_Send(send_bf.data(), count, MPI_LONG_LONG, 0, 3, MPI_COMM_WORLD);
        }
        MPI_Bcast(&min_butterflies, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
        if (min_butterflies == LLONG_MAX) {
            if (rank == 0) {
                cout << "\nIteration " << iter << ": No butterflies remaining\n";
            }
            break;
        }
        if (rank == 0) {
            sort(all_butterflies.begin(), all_butterflies.end(), 
                 [](const EdgeButterfly& a, const EdgeButterfly& b) {
                     return a.u == b.u ? a.v < b.v : a.u < b.u;
                 });
            long long total_butterflies = 0;
            for (const auto& eb : all_butterflies) {
                total_butterflies += eb.butterflies;
            }
            cout << "\nIteration " << iter << " Butterfly Counts:\n";
            cout << "----------------------------------------\n";
            cout << "| " << left << setw(6) << "U" << " | " << setw(6) << "V" << " | " << setw(12) << "Butterflies" << " | " << setw(10) << "Partition" << " |\n";
            cout << "----------------------------------------\n";
            for (const auto& eb : all_butterflies) {
                if (eb.butterflies > 0) {
                    int partition = 0;
                    for (int r = 0; r < size; ++r) {
                        if (find(local_graph.local_u.begin(), local_graph.local_u.end(), eb.u) != local_graph.local_u.end()) {
                            partition = rank;
                            break;
                        }
                    }
                    cout << "| " << left << setw(6) << eb.u << " | " << setw(6) << eb.v << " | " << setw(12) << eb.butterflies << " | " << setw(10) << partition << " |\n";
                }
            }
            cout << "----------------------------------------\n";
            cout << "  Total Butterflies: " << total_butterflies << "\n";
            cout << "  Minimum Butterfly Count: " << min_butterflies << "\n";
        }
        graph = peel_edges(graph, all_butterflies, min_butterflies, rank, size);
        MPI_Bcast(&graph.num_edges, 1, MPI_INT, 0, MPI_COMM_WORLD);
        u_offset.resize(graph.num_u + 1);
        u_edges.resize(graph.num_edges);
        if (rank == 0) {
            u_offset = graph.u_offset;
            u_edges = graph.u_edges;
        }
        MPI_Bcast(u_offset.data(), graph.num_u + 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(u_edges.data(), graph.num_edges, MPI_INT, 0, MPI_COMM_WORLD);
        if (rank != 0) {
            graph.u_offset = u_offset;
            graph.u_edges = u_edges;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    if (rank == 0) {
        double end_time = MPI_Wtime();
        cout << "\nFinal Summary:\n";
        cout << "  Remaining Edges: " << graph.num_edges << "\n";
        cout << "  Execution Time: " << fixed << setprecision(3) << (end_time - start_time) << " seconds\n";
        cout << "====================================\n";
    }
    MPI_Finalize();
    return 0;
}