#include <iostream>
#include "simulation.h"
#include <string.h>

using namespace std;

int main(int argc, char* argv[]) {

    Simulation sim;
    Rloop_equilibrium_model model;
    bool sandbox = false;
    //process command line arguments
    sim.set_infile(argv[1]);
    sim.set_outfile(argv[2]);
    for (int i=3; i<argc; i++) {
        if (!strcmp(argv[i], "--a")) {
            model.seta(atof(argv[i+1]));
            i++;
        }
        else if (!strcmp(argv[i], "--N")) {
            if (!strcmp(argv[i+1],"auto")){
                sim.set_auto_domain_size(true);
            }
            else{
                model.setN(atoi(argv[i+1]));
            }
            i++;
        }
        else if (!strcmp(argv[i], "--sigma")) {
            model.set_superhelicity(atof(argv[i+1]));
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
        else if (!strcmp(argv[i], "--invert")) {
            sim.complement_input();
            sim.reverse_input();
        }
        else if (!strcmp(argv[i], "--homopolymer")) {
            model.set_bp_energy_override(atof(argv[i+1]));
            i++;
        }
            //options specific to ensemble analyzer
        else if (!strcmp(argv[i], "--sandbox")) {
            sandbox = true;
        }
        else if (!strcmp(argv[i], "--bedfile")) {
            sim.set_bedfile(true);
        }
        else if (!strcmp(argv[i], "--unconstrained")) {
            model.set_unconstrained(true);
        }
        else if (!strcmp(argv[i], "--circular")) {
            sim.set_circular();
        }
        else if (!strcmp(argv[i], "--residuals")) {
            sim.set_residuals(true);
        }
        else if (!strcmp(argv[i], "--sensitivity")) {
            sim.set_power_threshold(atoi(argv[i+1]));
            i++;
        }
        else if (!strcmp(argv[i], "--top")) {
            sim.set_top(atoi(argv[i+1]));
            i++;
        }
        else if (!strcmp(argv[i], "--dump")) {
            sim.set_dump(true);
        }
        else if (!strcmp(argv[i], "--localaverageenergy")) {
            sim.set_average_g(true);
        }
        else{
            cout << "Unrecognized command line option: " << argv[i] << endl;
            exit(1);
        }
    }
    sim.add_model(model);
    if (sandbox){
        //sim.sandbox();
        sim.simulation_D();
        return 0;
    }
    sim.simulation_A();

    return 0;
}