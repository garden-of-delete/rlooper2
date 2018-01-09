#include <iostream>
#include "simulation.h"
#include <string.h>

using namespace std;

int main(int argc, char* argv[]) {

    Simulation sim;
    Rloop_equilibrium_model model;
    //process command line arguments
    sim.set_infile(argv[1]);
    sim.set_outfile(argv[2]);
    for (int i=3; i<argc; i++) {
        if (!strcmp(argv[i], "--sigma")) {
            model.set_superhelicity(atof(argv[i+1]));
            i++;
        }
        else if (!strcmp(argv[i], "--a")) {
            model.seta(atof(argv[i+1]));
            i++;
        }
        else if (!strcmp(argv[i], "--minlength")) {
            sim.set_minlength(atoi(argv[i+1]));
            i++;
        }
        else if (!strcmp(argv[i], "--reverse")) {
            sim.reverse_input();
        }
        else if (!strcmp(argv[i], "--complement")) {
            sim.complement_input();
        }
        else if (!strcmp(argv[i], "--sandbox")) {
            sim.add_model(model);
            sim.sandbox();
            return 0;
        }
        else if (!strcmp(argv[i], "--bedfile")) {
            sim.set_outfile2(argv[i+1]);
            i++;
        }
        else{
            cout << "Unrecognized command line option: " << argv[i] << endl;
            exit(1);
        }
    }
    sim.add_model(model);
    sim.simulation_A();

    return 0;
}