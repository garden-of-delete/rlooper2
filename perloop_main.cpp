//
// Created by Robert Stolz on 9/28/17.
//

#include <iostream>
#include "simulation.h"
#include <string.h>

using namespace std;

int main(int argc, char* argv[]) {

    //initialize new simulation
    Simulation sim;
    Rloop_equilibrium_model model;
    //print usage statement
    if (argc == 1){
        cout << "Usage: perloop sequence_file output_file superhelicity_lower_bound superhelicity_upper_bound" << endl;
        return 0;
    }
    //process command line arguments and set simulation parameters
    sim.set_infile(argv[1]);
    sim.set_outfile(argv[2]);
    float supercoiling_lower_bound = atof(argv[3]);
    float supercoiling_upper_bound = atof(argv[4]);
    for (int i=5; i<argc; i++) {
        if (!strcmp(argv[i], "--a")) {
            model.seta(atof(argv[i + 1]));
            i++;
        }
        else if (!strcmp(argv[i], "--minlength")) {
            sim.set_minlength(atoi(argv[i + 1]));
            i++;
        }
        else if (!strcmp(argv[i], "--reverse")){
            sim.reverse_input();
        }
        else if (!strcmp(argv[i], "--complement")){
            sim.complement_input();
        }
        else{
            cout << "Unrecognized command line option: " << argv[i] << endl;
            exit(1);
        }
    }

    //run simulation
    sim.add_model(model);
/*
    vector<vector<Structure>*> outer;
    vector<Structure>* inner = new vector<Structure>;
    inner->push_back(Structure());
    inner->push_back(Structure());
    inner->push_back(Structure());
    inner->push_back(Structure());
    inner->push_back(Structure());
    outer.push_back(inner);
    cout << outer.size() << ' ' << inner->size()<<endl;
    for (auto it = outer.begin(); it < outer.end(); ++it){
        delete *it;
    }
    outer.clear();
    cout << outer.size() << ' ' << inner->size()<<endl;
*/
        for (float superhelicity=supercoiling_lower_bound; superhelicity <= supercoiling_upper_bound+0.0001; superhelicity += 0.001){
            cout << superhelicity << endl;
            sim.simulation_B(superhelicity);
        }
    return 0;
}