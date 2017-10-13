//
// Created by Robert Stolz on 6/28/17.
//

#ifndef RLOOPER2_SIMULATION_H
#include "Rloop_equilibrium_model.h"
#import "exception_handling.h"
#import <sstream>

class Simulation{
private:
    std::vector<Model*> models;
    std::vector<Gene*> genes;
    ifstream infile;
    ofstream outfile;
    int minlength; //a minimum loop length to be applied simulation-wide.
    bool reverse_flag;
    bool complement_flag;
    //member functions

    /**
     *
     * @param gene
     */
    void write_wigfile(Gene& gene);

    /**
    * A test environmnet for debugging purposes
    */
    void sandbox();

public:
    Simulation();
    Simulation(int argc, char* argv[]);
    ~Simulation();
    //getters and setters
    void set_infile(string infilename);
    void set_outfile(string outfilename);
    void set_minlength(int Minlength);
    void reverse_input();
    void complement_input();
    std::vector<Model*> get_models();
    void add_model(Model& model);

    //simulation protocols
    /**
     * Ensemble analysis on a set of genes. Supports one model.
     */
    void simulation_A();

    /**
     * Computes P(R-Loop is on the sequence) for a given superhelicity level. Supports any number of models.
     */
    void simulation_B(float superhelicity);

    /**
     * Computes expected length for the ensemble at the given superhelicity value
     * @param superhelicity
     */
    void simulation_C(float superhelicity);
};

#define RLOOPER2_SIMULATION_H

#endif //RLOOPER2_SIMULATION_H
