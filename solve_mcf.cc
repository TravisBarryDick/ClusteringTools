#include <iostream>
#include <lemon/dimacs.h>
#include <fstream>
#include <lemon/smart_graph.h>
#include <lemon/capacity_scaling.h>
using namespace lemon;



int main(int argc, char *argv[]){
    //assert command line args
    if (argc < 3) {
        std::cerr << "Correct arguments: input_filename output_filename" << std::endl;
    }
    //read input
    std::ifstream f(argv[1]);
    SmartDigraph g;
    SmartDigraph::ArcMap<int> lower(g), capacity(g);
    SmartDigraph::ArcMap<float> cost(g);
    SmartDigraph::NodeMap<int> supply(g);
    readDimacsMin(f, g, lower, capacity, cost, supply);

    std::cerr << "Nodes: " << countNodes(g) << "\t Edges: " << countArcs(g) << std::endl;

    CapacityScaling<SmartDigraph> algo(g);
    algo.upperMap(capacity).costMap(cost).supplyMap(supply);
    
    //run teh algo
    algo.run();
    std::cerr << "Total cost: " << algo.totalCost() << std::endl;
    std::ofstream fout(argv[2]);

    //generate output
    for (SmartDigraph::ArcIt e(g) ; e != INVALID ; ++e) {
        if (algo.flow(e) != 0) {
            fout <<  g.id(g.source(e)) << " " << g.id(g.target(e)) << " " << algo.flow(e) << std::endl;
        }
    }
}
