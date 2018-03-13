//
// Created by Robert Stolz on 6/28/17.
//

#include "windower.h"

Windower::Windower(): min_window_size(2) {};

Windower::Windower(std::vector<char> &target_sequence){
    current_sequence = &target_sequence;
    current_start = target_sequence.begin();
    current_stop = current_start + min_window_size-1;
    circular = false;
}

int Windower::get_min_window_size(){
    return min_window_size;
}

void Windower::set_min_window_size(int size){
    if (size < 2){
        throw WindowerException();
    }
    min_window_size = size;
}

bool Windower::get_circular(){
    return circular;
}

void Windower::set_circular(bool value){
    circular = value;
}

void Windower::set_sequence(std::vector<char>& target_sequence){ //not working?
    current_sequence = &target_sequence;
    current_start = target_sequence.begin();
    current_stop = current_start + min_window_size - 2;
}

bool Windower::has_next_window(){
    //if the last window has been reached
    if (current_start == current_stop-min_window_size+1 && current_stop == current_sequence->end()-1){ //untested rs 7.11.17
        return false;
    }
    else
        return true;
}

void Windower::next_window_from_all_windows(std::vector<char>::iterator& start, std::vector<char>::iterator& stop){
    //if (!has_next_window()){ //safety check to make sure the next window exists
    //    throw WindowerException(); //throw exception?
    //}
    if (stop < current_sequence->end()-1){
        ++current_stop;
    }
    else{ //if (start < current_sequence->end()-min_window_size){
        ++current_start;
        current_stop = current_start + min_window_size-1;
    }
    start = current_start;
    stop = current_stop;
    return;
}

long int Windower::get_current_start_offset(){
    return current_start - current_sequence->begin();
}

long int Windower::get_current_stop_offset(){
    return current_stop - current_sequence->begin();
}

void Windower::reset_window(){
    current_start = current_sequence->begin();
    current_stop = current_sequence->begin()+min_window_size - 2;
}

void Windower::print_current_window(){
    std::vector<char>::iterator it;
    it = current_start;
    while (it < current_stop){
        std::cout << *it;
        ++it;
    }
    std::cout << '\n';
}