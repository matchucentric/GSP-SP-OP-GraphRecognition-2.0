//usage: ./dfs2 <input_graph>

#include <atomic>
#include <vector>
#include <iostream>
#include <fstream>
#include <stack>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <execution>

using std::pair;
using std::vector;

bool pairEqualUnordered(pair<int, int> p1, pair<int, int> p2) {
    // Small helper function to test equivalence of pairs, independant of order

    return p1.first == p2.first && p1.second == p2.second || p1.first == p2.second && p1.second == p2.first;
}

struct graphDataHolder{
    int n; //equal to |V| size of graph
    int e; //equal to |E| number of edges in the graph
    vector<int> parent; //vector containing data so that parent[3] gives the number of the vertex that is parent of 3 in the dfs tree
    vector<int> nDescendants; //vector containing number of descendants of a vertex including itself as a descendant.
    vector<int> dfsRank; //holds the number that represents in the order that the dfs accessed each vertex in the graph in.
    vector<int> adjList; //adjList[0] contains the vector of Edges that are connected to vertex 0. Edge contains target vertex and tree edge boolean.
    vector<int> adjAddress; //Holds address locations for each vector in adjlist
    vector<pair<int, int>> backEdges; //vector of back edges found in the dfs search
    vector<pair<int, int>> earVerticies; //lists the nodes of the back edge associated with the ear each node is in
    // pair.first must be source, pair.second must be sink. 
    // These label every vertex in the graph to a back edge
    vector<vector<int>> ear; //Holds complete list of edges for each ear.
};

graphDataHolder gD; //gD stands for Graph Data

bool isAncestor(int const & a, int const & b){
   return (gD.dfsRank[a] <= gD.dfsRank[b] && gD.dfsRank[b] < ( gD.dfsRank[a] + gD.nDescendants[a] ) );
}

bool lexiCompareVals(int const & p, int const & q, int const & x, int const & y){
    //return true if back-edge (q <-- p) is smaller than (y <-- x)
    //return false if (y <-- x) is smaller than (q <-- p)

    return x == -1 || (gD.dfsRank[q] < gD.dfsRank[y]) 
            || ( (gD.dfsRank[q] == gD.dfsRank[y]) && (gD.dfsRank[p] < gD.dfsRank[x]) && !isAncestor(p, x) ) 
            || ( (gD.dfsRank[q] == gD.dfsRank[y]) && isAncestor(x, p) );
}

bool lexiCompare(pair<int, int> p1, pair<int, int> p2) {
    // Wrapper function that takes 2 points for lexigraphical comparison of back edges
    // this version can be used in the std::sort
    return lexiCompareVals(p1.first, p1.second, p2.first, p2.second);
}

bool isTree(int const & a, int  const & b, graphDataHolder const & gD){
    // returns true iff a -> b or b -> a is a tree edge
    return (gD.parent[a] == b || gD.parent[b] == a);
}

namespace s=std;

void create_adjacency_list(int const &, vector<pair<int, int>> const &, vector<int> &, vector<int> &);
void create_edges(int const &, int const &, vector<pair<int, int>>&, s::ifstream&);
bool isAncestor(int  const &, int  const &);
//bool lexiCompare(vector<pair<int, int>>&, vector<pair<int, int>>&);
bool isTree(int const &, int const &, graphDataHolder const &);
void genCS(int const &);

int main(int argc, char* argv[]) {
    s::ios::sync_with_stdio(false);     
    
    // Pulling the first line defining number of nodes and edges from input files
    s::ifstream is(argv[1]);
    is >> gD.n >> gD.e;
    
    gD.dfsRank = s::vector<int>(gD.n, -1);
    gD.dfsRank.shrink_to_fit();
    gD.parent = s::vector<int>(gD.n, -1);
    gD.parent.shrink_to_fit();
    gD.nDescendants = s::vector<int>(gD.n, 1);
    gD.nDescendants.shrink_to_fit();
    gD.earVerticies = s::vector<pair<int, int>>(gD.n, s::pair<int, int>(-1, -1));
    gD.earVerticies.shrink_to_fit();
    gD.adjList.reserve(gD.e * 2);
    gD.adjAddress.reserve(gD.n+1);
    
    vector<pair<int, int>> edges;

    create_edges(gD.n, gD.e, edges, is);
    create_adjacency_list(gD.n, edges, gD.adjList, gD.adjAddress);
    gD.backEdges.reserve(gD.e - gD.n + 1);
    
    //measure the start time with a barrier that prevents the compiler from reordering
    std::atomic_thread_fence(std::memory_order_seq_cst);
    auto time_start = s::chrono::steady_clock::now();
    std::atomic_thread_fence(std::memory_order_seq_cst);
 
    genCS(0);
    
    //measure the end time with a barrier that prevents the compiler from reordering
    std::atomic_thread_fence(std::memory_order_seq_cst);
    auto time_end = s::chrono::steady_clock::now();
    std::atomic_thread_fence(std::memory_order_seq_cst);
    
    s::chrono::duration<double, s::nano> time_elapsed = s::chrono::duration_cast<s::chrono::nanoseconds>(time_end - time_start);
    s::cout << "time elapsed is: " << s::setprecision(4) << time_elapsed.count()/1000 << "μs\n";
    

    //Print adjacency list for debugging
    for(int i = 0; i<gD.n;i++){
        s::cout << "List " << i << " is: ";
        for(int j = gD.adjAddress[i]; j<gD.adjAddress[i+1]; j++){
            char c = (isTree(gD.adjList[j], i, gD)?'T':'B');
            s::cout << c << gD.adjList[j] << " ";
        }
        s::cout << '\n';
    }
    s::cout << '\n';

    //Print DFS ranks for testing
    s::cout << '\n';
    for (int i = 0; i < gD.dfsRank.size(); ++i) {
        s::cout << "Vertex " << i << " | DFS " << gD.dfsRank[i] << " | Parent " << gD.parent[i] << " | nDescendants " << gD.nDescendants[i] << " | ear " << gD.earVerticies[i].first << ',' << gD.earVerticies[i].second << "\n";
    }

    //Print Ears for testing
    s::cout << '\n';
    for (int i = 0; i < gD.ear.size(); i++) {
        s::cout << "Ear " << i << " | Nodes: ";
        vector<int> ear = gD.ear[i];
        for (int j = 0; j < ear.size(); j++) {
            s::cout << ear[j] << " ";   
        }
        s::cout << "\n"; 
    }

    // return 0;
}

void genCS(int const & starting_vertex) {

    // The stack holds vertex numbers. they are popped off when all edges out from
    // the vertex have been visited and all decendants have as well
    std::vector<int> the_stack;

    // DFS keeps track of visited verticies while this vector holds edges that have been visited
    // so they cannot be again
    std::vector<pair<int, int>> visited_edges;

    the_stack.reserve(gD.n);
    the_stack.reserve(gD.e * gD.e);
    the_stack.push_back(starting_vertex);

    int dfsNumber = 1;
    gD.dfsRank[starting_vertex] = dfsNumber++;

    while (!the_stack.empty()) {

        // current node being visited by dfs
        int topOfStack = the_stack.back();
        // descend becomes true when a new, unvisited node is found. it is added
        // to the stack and immediately starts being visited in next iteration
        bool descend = false;

        for (int i = gD.adjAddress[topOfStack]; i < gD.adjAddress[topOfStack+1]; i++) {
            // Current edge being viewed is topOfStack -> w
            int w = gD.adjList[i];

            // repeated is empty if the edge has not been visited. uses
            // custom comparison function
            auto repeated = std::find_if(visited_edges.begin(), visited_edges.end(), 
            [&](pair<int, int> p) { return (pairEqualUnordered(p, std::make_pair(topOfStack, w))); });

            // skips the edge if its been visited
            if (repeated != visited_edges.end()) {
                continue;
            }

            // Add edge to the list of visited ones
            visited_edges.push_back(std::make_pair(topOfStack, w));


            if (gD.dfsRank[w] == -1) {
                // New tree edge detected

                // Set the new nodes parent
                gD.parent[w] = topOfStack;
                // Set the new nodes dfs number
                gD.dfsRank[w] = dfsNumber++;
                // add new node to the stack
                the_stack.push_back(w);
                // set flag for new iteration and proceed there now
                descend = true;
                break;

            }else if(gD.dfsRank[w] < gD.dfsRank[topOfStack] && w != gD.parent[topOfStack]){
                //back edge detected

                // Add the back edge to the list for making ears later
                gD.backEdges.push_back(std::make_pair(topOfStack, w));

                // If the node is currently not part of an ear or the new back edge is
                // lexigraphically before the previous one associated with this node
                if (gD.earVerticies[topOfStack].first == -1 || !lexiCompare(gD.earVerticies[topOfStack], 
                std::make_pair(w, topOfStack))) {
                    
                    // Set the back edge for this node to this edge
                    gD.earVerticies[topOfStack].first = topOfStack;
                    gD.earVerticies[topOfStack].second = w;

                }

                
            }
        }

        // descend to new ndoe and ignore rest of iteration
        if (descend) {
            continue;
        }
        
        // Recursively adding descendants to parents once all tree edges have been found
        if(topOfStack != starting_vertex){
            gD.nDescendants[gD.parent[topOfStack]] += gD.nDescendants[topOfStack];
        }
        /*
        dfsRank[topOfStack] > 2; this is an easier way of writing
        topOfStack != starting_vertex && parent[topOfStack] != starting_vertex
        
        both of those are necessary because:
        1. if I only check topOfStack != starting_vertex, then when the vertex with the root r as parent is backed into as the searVerticiesch
            is reversing up the tree, the if condition will succeed and gD.earVerticies[gD.parent[topOfStack]] will run which produces an error
            because the result is -1 and the function will use that parameter to find gD.dfsRank[-1] which doesn't exist.
        2. if I only check parent[topOfStack] != starting_vertex, then when the root r is backed into as the searVerticiesch reverses up the tree
            the if condition will succeed and gD.earVerticies[gD.parent[topOfStack]] will run which produces an error because gD.earVerticies[-1] doesn't
            exist.
        */
        int current = topOfStack;

        // Recursively travelling up the tree to replace the nodes corrosponding back edges
        // to the one associated with this node if it is lexigraphically smaller.
        while (gD.dfsRank[current] > 2) {
 
            if(lexiCompare(gD.earVerticies[current], gD.earVerticies[gD.parent[current]])){
                    gD.earVerticies[gD.parent[current]].first = gD.earVerticies[current].first; 
                    gD.earVerticies[gD.parent[current]].second = gD.earVerticies[current].second;
            } else {
                break;
            }

            current = gD.parent[current];
        }

        the_stack.pop_back();
    }

    // sort the back edges with custom comparison function, so they are in order of their ear number
    sort(gD.backEdges.begin(), gD.backEdges.end(), lexiCompare);

    // adding edges to each ear, starting with the back edge until it reaches a node that does
    // not corrospond with the back edge
    for (auto const& edge : gD.backEdges) {

        // creating the ear vector to be added to master list
        vector<int> ear;
        ear.reserve(gD.n + 1);

        // adding the two back edge nodes to start
        ear.push_back(edge.second);
        ear.push_back(edge.first);
        // start with higher dfs number node and work back up the tree
        int current = edge.first;

        while (true) {
            
            // retreiving back edge for current node to test against back edge for the current node
            pair<int, int> backEdge = gD.earVerticies[current];
            if (backEdge.first == edge.first && backEdge.second == edge.second) {

                // push parent on (as ear must contain the next node if the current vertex belongs to this ear)
                current = gD.parent[current];
                ear.push_back(current);

            } else {
                // end of ear
                break;
            }
        }
        // constructing ears back edge first was simple, but they should appear with the back edge last
        std::reverse(ear.begin(), ear.end());
        ear.shrink_to_fit();
        // add ear to master list
        gD.ear.push_back(ear);
    }
}

void create_adjacency_list(int const & n, vector<pair<int, int>> const & edges, vector<int>& adj, vector<int> & adjAdd) {

    for (int i = 0; i < n + 1; i++) {
        adjAdd.emplace_back(0);
    }

    for (int i = 0; i < edges.size() * 2; i++) {
        adj.emplace_back(0);
    }

    vector<int> temp(n+1, 0);
    for(int i = 0; i < edges.size(); i++){
        adjAdd[edges[i].first+1]++;
        adjAdd[edges[i].second+1]++;
    }
    
    for(int i = 1; i<=n; i++){
        adjAdd[i] += adjAdd[i-1];
        temp[i] = adjAdd[i];
    }
    
    for (int i = 0; i<edges.size(); i++) {
        int x = edges[i].first;
        int y = edges[i].second;
        adj[temp[x]++] = y;
        adj[temp[y]++] = x;
    }
}

void create_edges(int const & num_v, int const & num_e, vector<pair<int, int>>& edges, s::ifstream& istream) {
    edges.reserve(num_e);
    for (int i = 0; i < num_e; ++i) {
        int x, y;
        istream >> x >> y;
        edges.emplace_back( (x>y?y:x), (x>y?x:y));
    }
}