//
// Created by Robert Stolz on 6/28/17.
//

#ifndef RLOOPER2_SIMULATION_H
#include "Rloop_equilibrium_model.h"
#import "exception_handling.h"
#import <sstream>

class Simulation{
private:
    std::vector<Model> models;
    std::vector<Gene*> genes;
    ifstream infile;
    ofstream outfile;
    //member functions
    void write_wigfile(Gene& gene);
    //simulation protocols

    /**
     * Ensemble analysis on a set of genes. Further description needed.
     */
    void simulation_A();

    /**
     * Computes P(R-Loop is on the sequence at equilibrium for a given supercoiling level)
     */
    void simulation_B(float superhelicity);

    /**
    * A test environmnet for debugging purposes
    */
    void sandbox();

public:
    Simulation(int argc, char* argv[]);
    ~Simulation();
    //getters and setters
    void set_infile(string infilename);
    void set_outfile(string outfilename);

    /**
     * Tells the software suites what simulations to run in what order. Functions as the main method.
     */
    void run_simulations(); //main method for the software suite.

};

#define RLOOPER2_SIMULATION_H

#endif //RLOOPER2_SIMULATION_H
