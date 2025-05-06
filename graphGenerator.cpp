#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>

using namespace std;

int main() {
    int num_u = 128;
    int num_v = 707;
    int total_edges = 708;
    vector<vector<int>> adj_list(num_u);
    adj_list[0].resize(122);
    for (int i = 0; i < 122; ++i) {
        adj_list[0][i] = i * 5 % num_v;
    }
    int remaining_edges = total_edges - 122;
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> edge_dist(1, 11);
    vector<int> u_indices(num_u - 1);
    for (int i = 0; i < num_u - 1; ++i) {
        u_indices[i] = i + 1;
    }
    shuffle(u_indices.begin(), u_indices.end(), gen);
    int edges_assigned = 0;
    for (int i = 0; i < num_u - 1 && edges_assigned < remaining_edges; ++i) {
        int u = u_indices[i];
        int max_edges = min(edge_dist(gen), remaining_edges - edges_assigned);
        adj_list[u].resize(max_edges);
        uniform_int_distribution<> v_dist(0, num_v - 1);
        for (int j = 0; j < max_edges; ++j) {
            adj_list[u][j] = v_dist(gen);
        }
        sort(adj_list[u].begin(), adj_list[u].end());
        adj_list[u].erase(unique(adj_list[u].begin(), adj_list[u].end()), adj_list[u].end());
        edges_assigned += adj_list[u].size();
    }
    while (edges_assigned < remaining_edges) {
        int u = u_indices[0];
        adj_list[u].push_back(uniform_int_distribution<>(0, num_v - 1)(gen));
        sort(adj_list[u].begin(), adj_list[u].end());
        adj_list[u].erase(unique(adj_list[u].begin(), adj_list[u].end()), adj_list[u].end());
        edges_assigned = 0;
        for (int i = 1; i < num_u; ++i) {
            edges_assigned += adj_list[i].size();
        }
    }
    ofstream out("graph.txt");
    out << "AdjacencyGraph\n";
    out << num_u << "\n";
    out << total_edges << "\n";
    for (int u = 0; u < num_u; ++u) {
        for (size_t i = 0; i < adj_list[u].size(); ++i) {
            out << adj_list[u][i];
            if (i < adj_list[u].size() - 1) {
                out << " ";
            }
        }
        out << "\n";
    }
    out.close();
    return 0;
}