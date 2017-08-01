//
// Created by Robert Stolz on 6/28/17.
//

#ifndef RLOOPER2_WINDOWER_H
#define RLOOPER2_WINDOWER_H

#include <vector>

class Windower {
private:
    int min_window_size;
    std::vector<char>* current_sequence;
    std::vector<char>::iterator current_start, current_stop;
public:
    //constructors
    Windower();
    Windower(std::vector<char> &sequence);

    void set_sequence(std::vector<char>& sequence);

    /**
     * Peeks the next window of the ?current windowing scheme? and returns the result as a boolean
     * @return  if true, the next window exists, if false, it does not.
     */
    bool has_next_window();

    void next_window_from_all_windows(std::vector<char>::iterator& start, std::vector<char>::iterator& stop);

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
};


#endif //RLOOPER2_WINDOWER_H
