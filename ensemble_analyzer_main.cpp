#include <iostream>
#include "simulation.h"

using namespace std;

int main(int argc, char* argv[]) {

    Simulation sim;
    Rloop_equilibrium_model model;
    //process command line arguments
    sim.set_infile(argv[1]);
    sim.set_outfile(argv[2]);
    for (int i=3; i<argc; i++){
        if (!strcmp(argv[i],"--sigma")){
            model.setSigma(atof(argv[i+1]));
            i++;
        }
        else if (!strcmp(argv[i],"--a")){
            model.seta(atof(argv[i+1]));
            i++;
        }
        else if (!strcmp(argv[i],"--minlength")){
            sim.set_minlength(atoi(argv[i+1]));
            i++;
        }

    }
    sim.simulation_A();

    return 0;
}