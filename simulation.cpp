//
// Created by Robert Stolz on 6/28/17.
//
#import "simulation.h"

Simulation::Simulation(){

}

Simulation::~Simulation(){
    for(std::vector<Gene*>::iterator it = genes.begin(); std::distance(genes.begin(),it) < genes.size(); ++it){
        delete *it; //need to test this destructor
    }
}

void Simulation::set_infile(string infilename){
    infile.open(infilename, ios::in);
}

void Simulation::set_outfile(string outfilename){
    outfile.open(outfilename, ios::out);
}

void Simulation::write_wigfile(Gene& gene){

    //compute the r-loop involvement probability for each base (will probably be moved out of this func later)
    std::vector<double> bp_probabilities(gene.get_length(),0.0);
    //for each structure in the gene
    for (std::vector<Structure>::iterator it = gene.getStructures()[0]->begin();
         it < gene.getStructures()[0]->end(); ++it){
        //for each base in the structure
        for(long int i=it->position.start_pos-gene.getPosition().start_pos;
            i <it->position.end_pos-gene.getPosition().start_pos; i++){
            bp_probabilities[i] += it->probability;
        }
    }
    //if strand is -, reverse bp_probabilities
    if (gene.getPosition().strand == "-"){
        std::reverse(bp_probabilities.begin(),bp_probabilities.end());
    }
    //open filestream
    std::stringstream ss;
    string wigfile_name = gene.getHeader().c_str();
    //wigfile_name = wigfile_name + ".wig";
    wigfile_name = gene.getHeader() + ".wig";
    ofstream wigfile(wigfile_name,ios::out);
    //compose .wig header
    string name = gene.getName();
    //adjust browser position
    ss << "browser position " << gene.getPosition().chromosome << ':' << gene.getPosition().start_pos << '-' <<
       gene.getPosition().end_pos << endl;
    ss << "track type=wiggle_0 name=\"" << name << "\" visibility=full autoscale=off color=50,150,255 priority=10" << endl;
    ss << "fixedStep chrom=" << gene.getPosition().chromosome << " start=" << gene.getPosition().start_pos << " step=1" << endl;
    for (int i=0; i < bp_probabilities.size(); i++){
        ss << bp_probabilities[i] << endl;
    }
    //write stringstream to file
    wigfile << ss.rdbuf();
}

void Simulation::simulation_A(){ //some of this code might be migrated into new objects and functions in the future
    //initialize variables
    if (!infile.is_open()){
        //throw exception
    }
    bool eof = false;
    Rloop_equilibrium_model modelA;

    //do while !eof
    while(eof == false){
        //allocate new gene
        Gene* this_gene = new Gene();
        //read gene
        eof = this_gene->read_gene(infile);
        //compute structures using models
        this_gene->complement_sequence();
        this_gene->compute_structures(modelA);

        //ensemble analysis, free energies and boltzmann factors have already been computed in compute_structures
        //compute partition function
        long double partition_function = 0;
        long double sanity_check = 0;
        for (vector<Structure>::iterator it = this_gene->getStructures()[0]->begin();
             it < this_gene->getStructures()[0]->end(); ++it){
               partition_function += it->boltzmann_factor;
        }
        partition_function += modelA.ground_state_factor();
        //compute boltzmann weights and store in the structures
        for (vector<Structure>::iterator it = this_gene->getStructures()[0]->begin();
             it < this_gene->getStructures()[0]->end(); ++it){
            it->probability = it->boltzmann_factor/partition_function;
            sanity_check += it->boltzmann_factor/partition_function;
        }
        //compute p(base-pair i is in an R-Loop structure) and write to file
        write_wigfile(*this_gene);

        //unload the sequence data from the gene to save memory
        this_gene->unload(); //untested
        //store the gene in the genes vector
        genes.push_back(this_gene);
    }

}

void Simulation::run_simulations(){ //main method
    simulation_A();
}
