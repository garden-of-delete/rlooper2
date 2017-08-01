//
// Created by Robert Stolz on 6/28/17.
//

#ifndef RLOOPER2_SIMULATION_H
#include "rloop_equilibrium_model.h"
#import "exception_handling.h"
#import <sstream>

class Simulation{
private:
    std::vector<Model> models;
    std::vector<Gene*> genes;

    //member functions
    void write_wigfile(Gene& gene);
    //simulation protocols
    /**
     * Simulation protocol and reasoning goes here.
     */
    void simulation_A();


public:
    Simulation();
    ~Simulation();

    /**
     * Tells the software suites what simulations to run in what order. Functions as the main method.
     */
    void run_simulations(); //main method for the software suite.

};

#define RLOOPER2_SIMULATION_H

#endif //RLOOPER2_SIMULATION_H
