#include <iostream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/degree_centrality.hpp>
// #include <boost/graph/closeness_centrality.hpp>
// #include <boost/graph/betweenness_centrality.hpp>

using namespace std;

struct Atom { // can be read from the gro file
    string parent_code;
    string atom_code;
    double x_coord, y_coord, z_coord;
    double x_vel, y_vel, z_vel;
};

int main() {
    // need to replace this with the actual bond data
    vector<pair<int, int>> bonds = {{1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, {7, 8}, {8, 9}, {9, 10}};

    using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, Atom>;
    Graph g(bonds.begin(), bonds.end(), 10);

    // Degree Centrality
    vector<double> degree_centrality(boost::num_vertices(g));
    boost::degree_centrality(g, boost::make_iterator_property_map(degree_centrality.begin(), boost::get(boost::vertex_index, g), degree_centrality[0]));

    // Closeness Centrality
    // vector<double> closeness_centrality(boost::num_vertices(g));
    // boost::closeness_centrality(g, boost::make_iterator_property_map(closeness_centrality.begin(), boost::get(boost::vertex_index, g), closeness_centrality[0]));

    // Betweenness Centrality
    // vector<double> betweenness_centrality(boost::num_vertices(g));
    // boost::betweenness_centrality(g, boost::make_iterator_property_map(betweenness_centrality.begin(), boost::get(boost::vertex_index, g), betweenness_centrality[0]));

    // Results
    for (size_t i = 0; i < degree_centrality.size(); ++i) {
        cout << "Node " << i + 1
             << ": Degree Centrality = " << degree_centrality[i]
            //  << ", Closeness Centrality = " << closeness_centrality[i]
            //  << ", Betweenness Centrality = " << betweenness_centrality[i]
             << endl;
    }

    return 0;
}
