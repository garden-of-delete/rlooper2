//
// Created by Robert Stolz on 6/28/17.
//
#import "simulation.h"

void initialize_models(char* model){

}

Simulation::Simulation(){
    //default constructor
    minlength = 20; //default minlength
    reverse_flag = false;
    complement_flag = false;
}

Simulation::~Simulation(){
    for(std::vector<Gene*>::iterator it = genes.begin(); std::distance(genes.begin(),it) < genes.size(); ++it){
        delete *it; //need to test this destructor
    }
    outfile.close();
    outfile2.close();
}

void Simulation::set_infile(string infilename){
    infile.open(infilename, ios::in);
}

void Simulation::set_outfile(string outfilename){
    outfile.open(outfilename, ios::out);
}

void Simulation::set_outfile2(string outfilename){
    outfile2.open(outfilename, ios::out);
}

void Simulation::set_minlength(int Minlength){
    minlength = Minlength;
}

void Simulation::reverse_input(){
    reverse_flag = true;
}

void Simulation::complement_input(){
    complement_flag = true;
}

std::vector<Model*> Simulation::get_models(){
    return models;
}

void Simulation::add_model(Model& model){
    models.push_back(&model);
}

void Simulation::compute_signal_bpprobs(Gene &gene, vector<double> *&signal){
    signal = new vector<double>(gene.get_length(), 0.0);
    //compute the r-loop involvement probability for each base (will probably be moved out of this func later)
    //for each structure in the gene
    for (std::vector<Structure>::iterator it = gene.getStructures()[0]->begin();
         it < gene.getStructures()[0]->end(); ++it) {
        //for each base in the structure
        for (long int i = it->position.start_pos - gene.getPosition().start_pos;
             i < it->position.end_pos - gene.getPosition().start_pos; i++) {
            (*signal)[i] += it->probability;
        }
    }
    //if strand is -, reverse bp_probabilities
    if (gene.getPosition().strand == "-") {
        std::reverse(signal->begin(), signal->end());
    }
}

void Simulation::call_peaks_threshold(Gene& gene, vector<double>& signal, vector<Loci>& peaks){
    int power_threshold = 6; //needs to be made a class variable
    double minimum = 1;
    bool in_peak = false;
    long peak_start=0, peak_end=0;
    double magnitude = 0;
    Structure* temp;
    for (int i=0; i < signal.size(); i++){
        //determine lowest value in the signal
        if (signal[i] < minimum){
            minimum = signal[i];
        }
    }
    for (int i=0; i < signal.size(); i++){
        if (signal[i] > minimum*pow(10,power_threshold)){ //the signal is significant
            if (!in_peak){
                in_peak = true;
                peak_start = gene.getPosition().start_pos + i;
            }
        }
        else{ //the signal is not significant
            if (in_peak){
                in_peak = false;
                peak_end = gene.getPosition().start_pos + i;
                peaks.emplace_back(Loci(gene.getPosition().chromosome,gene.getPosition().strand, peak_start, peak_end)); //chromosome, strand, start_pos, end_pos
            }
        }
    }
}

void Simulation::write_wigfile(Gene* gene, std::vector<double>* signal){
    //open stringstream
    std::stringstream ss;
    string wigfile_name = gene->getHeader().c_str();
    //compose .wig header
    string name = gene->getName();
    //adjust browser position
    ss << "browser position " << gene->getPosition().chromosome << ':' << gene->getPosition().start_pos << '-' <<
       gene->getPosition().end_pos << endl;
    ss << "track type=wiggle_0 name=\"" << name << "\" visibility=full autoscale=off color=50,150,255 priority=10"
       << endl;
    ss << "fixedStep chrom=" << gene->getPosition().chromosome << " start=" << gene->getPosition().start_pos << " step=1"
       << endl;
    for (int i = 0; i < signal->size(); i++) {
        ss << (*signal)[i] << endl;
    }
    //write stringstream to file
    outfile << ss.rdbuf();
}

void Simulation::write_bedfile(Gene* gene, vector<Loci>& peaks){
    //write bedfile
    stringstream ss;
    string strand_name;
    int start_pos=0, end_pos=0;

    if (gene->getPosition().strand == "+"){
        strand_name = "POS";
    }
    else {
        strand_name = "NEG";
    }
    ss << "browser position " << gene->getPosition().chromosome << ':' << gene->getPosition().start_pos << '-' <<
       gene->getPosition().end_pos << endl;
    ss << "track name=rLooper description=\""<< gene->getName()<<"\" useScore=1" << endl;
    //print BED header here
    //print the peaks in BED format
    for (int i=0; i < peaks.size(); i++){
        ss << peaks[i].chromosome << '\t' << (peaks)[i].start_pos << '\t' << peaks[i].end_pos
           << '\t' << strand_name << i << '\t' << '0' << '\t' << peaks[i].strand << endl;
    }
    //write stringstream to file
    outfile2 << ss.rdbuf();
}

void Simulation::simulation_A(){ //some of this code might be migrated into new objects and functions in the future
    //initialize variables
    if (!infile.is_open()){
        throw UnexpectedClosedFileException("Simulation::simulation_A");
    }
    bool eof = false;
    if (models.size() < 1){
        //throw exception
    }
    //do while !eof
    while(eof == false){
        //allocate new gene
        Gene* this_gene = new Gene();
        this_gene->windower.set_min_window_size(minlength);
        //read gene
        eof = this_gene->read_gene(infile);
        cout << "processing gene: " << this_gene->getName() << "...";
        //compute structures using models
        if (complement_flag)
            this_gene->complement_sequence();
        if (reverse_flag)
            this_gene->invert_sequence();
        this_gene->compute_structures(*models[0]);

        //ensemble analysis, free energies and boltzmann factors have already been computed in compute_structures
        //compute partition function
        long double partition_function = 0;
        long double sanity_check = 0;
        for (vector<Structure>::iterator it = this_gene->getStructures()[0]->begin();
             it < this_gene->getStructures()[0]->end(); ++it){
               partition_function += it->boltzmann_factor;
        }
        partition_function += models[0]->ground_state_factor();
        //compute boltzmann weights and store in the structures
        for (vector<Structure>::iterator it = this_gene->getStructures()[0]->begin();
             it < this_gene->getStructures()[0]->end(); ++it){
            it->probability = it->boltzmann_factor/partition_function;
            sanity_check += it->boltzmann_factor/partition_function;
        }
        sanity_check += models[0]->ground_state_factor()/partition_function;
        double test = fabs(1-sanity_check);
        if (fabs(1-sanity_check) > .00001){
            throw SimulationException("Probsum != 1"); //this throw is uncaught
        }
        vector<double>* signal = NULL;
        vector<Loci> peaks;
        compute_signal_bpprobs(*this_gene, signal);
        write_wigfile(this_gene,signal);
        //call peaks
        if (outfile2.is_open()){
            call_peaks_threshold(*this_gene,*signal,peaks); //possible null pointer exception generated here
            //write to bedfile
            write_bedfile(this_gene,peaks);
        }
        cout << "complete!" << endl;
        delete signal;
        //clear_sequence the sequence data from the gene to save memory
        this_gene->clear_sequence();
        this_gene->clear_structures();
        //store the gene in the genes vector
        genes.push_back(this_gene);
    }
}

//computes P(R-Loop) for the provided supercoiling value
void Simulation::simulation_B(float superhelicity){
    if (!infile.is_open()){
        throw UnexpectedClosedFileException("Simulation::simulation_B");
    }
    if (models.size() < 1){
        //throw exception
    }
    float p_rloop = 0;
    Gene* this_gene;
    if (!genes.size()){
        this_gene = new Gene();
        this_gene->read_gene(infile);
        this_gene->windower.set_min_window_size(minlength);
        this_gene->complement_sequence();
        //this_gene->invert_sequence();
        genes.push_back(this_gene);
    }
    else{
        this_gene = genes[0];
    }
    models[0]->set_superhelicity(superhelicity); //set the superhelicity in the model to the provided value
    this_gene->clear_structures(); //saves memory
    this_gene->compute_structures(*(models[0]));
    //determine P(ground state)
    long double partition_function = 0;
    long double ground_state_factor = 0;
    int index  = this_gene->getStructures().size();
    int count = 0;
    for (vector<Structure>::iterator it = this_gene->getStructures()[index-1]->begin();
         it < this_gene->getStructures()[index-1]->end(); ++it){
        partition_function += it->boltzmann_factor;
        count++;
    }
    ground_state_factor = models[0]->ground_state_factor();
    partition_function += ground_state_factor;
    //determine P(R-Loop) as 1-P(ground state)
    p_rloop = 1 - (ground_state_factor/partition_function);
    //display result
    outfile << superhelicity << ' ' << p_rloop << endl;
}

void Simulation::simulation_C(float superhelicity){
    if (!infile.is_open()){
        throw UnexpectedClosedFileException("Simulation::simulation_C");
    }
    if (models.size() < 1){
        //throw exception
    }
    Gene* this_gene;
    if (!genes.size()){
        this_gene = new Gene();
        this_gene->read_gene(infile);
        this_gene->windower.set_min_window_size(minlength);
        this_gene->complement_sequence();
        //this_gene->invert_sequence();
        genes.push_back(this_gene);
    }
    else{
        this_gene = genes[0];
    }
    models[0]->set_superhelicity(superhelicity); //set the superhelicity in the model to the provided value
    this_gene->compute_structures(*models[0]);
    //determine P(ground state)
    long double partition_function = 0;
    long double ground_state_factor = 0;
    int index  = this_gene->getStructures().size();
    int count = 0;
    for (vector<Structure>::iterator it = this_gene->getStructures()[index-1]->begin();
         it < this_gene->getStructures()[index-1]->end(); ++it){
        partition_function += it->boltzmann_factor;
        count++;
    }
    ground_state_factor = models[0]->ground_state_factor();
    partition_function += ground_state_factor;
    //determine expected length at the given superhelicity value
    double expected_length = 0;
    for (vector<Structure>::iterator it = this_gene->getStructures()[index-1]->begin();
         it < this_gene->getStructures()[index-1]->end(); ++it){
        expected_length += (it->boltzmann_factor/partition_function)*it->position.get_length();
    }
    //display result
    outfile << superhelicity << ' ' << expected_length << endl;
}

void Simulation::sandbox(){ //test/debug environment
    if (!infile.is_open()){
        throw UnexpectedClosedFileException("Simulation::sandbox");
    }
    Gene geneA;
    geneA.read_gene(infile); //sequence being read in is not used for anything
    //craig's simulation
    double lower_bound = -0.5;
    double upper_bound = 1.0;
    double supercoiling = 0.0;
    double last_supercoiling = 0.0;
    vector<double> x;
    vector<double> y;
    double step_size = .01;
    double tolerance = 0.05;
    char last_direction = 'n';

    Rloop_equilibrium_model modelA;
    //modelA.setMinimum_loop_length(minlength); //not functional, needs to be removed
    geneA.windower.set_min_window_size(minlength);
    //for each base pairing energy
    for (double bp_energy = lower_bound; bp_energy <= upper_bound; bp_energy += 0.1){
        //for each level of supercoiling
        supercoiling = last_supercoiling;
        while(true) {
            cout << "For bp_energy: " << bp_energy << ", and supercoiling: " << supercoiling << endl;
            //set supercoiling
            modelA.set_superhelicity(supercoiling);
            modelA.set_bp_energy_override(bp_energy);
            //dump previous ensemble of structures
            geneA.clear_structures();
            //compute ensemble
            geneA.compute_structures(modelA);
            //sum the probabilities of all R-loop structures
            double partition_function = 0;
            for (vector<Structure>::iterator it = geneA.getStructures()[0]->begin(); it != geneA.getStructures()[0]->end(); ++it){
                partition_function += it.base()->boltzmann_factor;
            }
            //compare to the partition function
            // if the P(R-loop) is > 50% + tol
            if (partition_function/(partition_function+modelA.ground_state_factor()) < .5){
                //adjust supercoiling down
                if (last_direction == 'u'){
                    //overshot the target and outside the error tolerance
                    step_size /= 0.5; //cut step size in half
                }
                supercoiling -= step_size;
                last_direction = 'd';
            }
                //if the {(R-loop) is < 50% - tol
            /*else if (partition_function/(partition_function+modelA.ground_state_factor()) < .5-tolerance){
                //adjust supercoiling up
                if (last_direction == 'd'){
                    //overshot the target and outside the error tolerance
                    step_size /= 0.5; //cut step size in half
                }
                supercoiling += step_size;
                last_direction = 'u';
            }*/
            else { //hit target within tolerance
                cout << "complete" << endl;
                x.push_back(bp_energy);
                y.push_back(supercoiling);
                last_supercoiling = supercoiling; //improves search efficiency at the next bp_energy
                break;
            }
        }
    }
    //write result
    stringstream ss;
    for (int i=0; i < x.size(); i++){
        ss << x[i] << ' ' << y[i] << endl;
    }
    outfile << ss.rdbuf();
    outfile.close();
}