#include <iostream>
#include "simulation.h"

using namespace std;

int main(int argc, char* argv[]) {

    Simulation sim(argc, argv);
    sim.run_simulations();
    return 0;
}