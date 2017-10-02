//
// Created by Robert Stolz on 9/28/17.
//

#include <iostream>
#include "simulation.h"

using namespace std;

int main(int argc, char* argv[]) {

    //initialize new simulation
    Simulation sim;
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

    //run simulation
        for (float superhelicity=supercoiling_lower_bound; superhelicity <= supercoiling_upper_bound+0.0001; superhelicity += 0.01){
            sim.simulation_B(superhelicity);
        }
    return 0;
}