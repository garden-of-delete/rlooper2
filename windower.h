//
// Created by Robert Stolz on 6/28/17.
//

#ifndef RLOOPER2_WINDOWER_H
#define RLOOPER2_WINDOWER_H

#include <vector>
#include <iostream>
#include "exception_handling.h"

class Windower {
private:
    long int min_window_size, current_window_size, sequence_size, sliding_window_size; //must be at least 2
    std::vector<char>* current_sequence;
    std::vector<char>::iterator current_start, current_stop;
    bool is_circular;
public:
    //constructors
    Windower();
    Windower(std::vector<char> &sequence);

    int get_min_window_size();

    void set_min_window_size(int size);

    void set_circular(bool Is_circular);

    void set_sequence(std::vector<char>& sequence);

    /**
     * Peeks the next window of the ?current windowing scheme? and returns the result as a boolean
     * @return  if true, the next window exists, if false, it does not.
     */
    bool has_next_window();

    bool has_next_window_circular();

    void next_window_from_all_windows(std::vector<char>::iterator& start, std::vector<char>::iterator& stop);

    void next_window_from_all_windows_circular(std::vector<char>::iterator& start, std::vector<char>::iterator& stop);

    void next_sliding_window(std::vector<char>::iterator& start, std::vector<char>::iterator& stop);

    void set_sliding_window_size(int size);

    void grow_window_right(std::vector<char>::iterator& start, std::vector<char>::iterator& stop);

    void grow_window_left(std::vector<char>::iterator& start, std::vector<char>::iterator& stop);

    void shrink_window_right(std::vector<char>::iterator& start, std::vector<char>::iterator& stop);

    void shrink_window_left(std::vector<char>::iterator& start, std::vector<char>::iterator& stop);

    bool next_sliding_window_length_n(std::vector<char>::iterator start, std::vector<char>::iterator stop);

    /**
     * Returns the offset of the current_start from the position where current_sequence begins
     * @return  a long int representing the offset
     */
    long int get_current_start_offset();

    /**
     * Returns the offset of the current_end from the position where current_sequence begins
     * @return  a long int representing the offset
     */
    long int get_current_stop_offset();

    /**
     * resets the current_start and current_stop iterators to their initial position on current_sequence.
     */
    void reset_window();

    /**
     * A function primarily for debugging. Allows the windower object to print the contents of the current window to stdout.
     *
     */
    void print_current_window();

};


#endif //RLOOPER2_WINDOWER_H
