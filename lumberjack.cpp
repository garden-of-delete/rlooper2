//
// Created by Robert Stolz on 8/8/17.
//

#include "lumberjack.h"

Lumberjack::Lumberjack(int _logging_level){
    logging_level = _logging_level;
    stringstream ss;
    ss << get_time() << "_log.txt";
    logfile.open(ss.str(), ios::out);
}

Lumberjack::~Lumberjack(){
    if (logfile.is_open()) {
        logfile.close();
    }
}

string Lumberjack::get_time(){
    auto t = time(nullptr);
    auto tm = *localtime(&t);

    std::ostringstream oss;
    oss << std::put_time(&tm, "%d-%m-%Y %H-%M-%S");
    string str = oss.str();
    return str;
}

void Lumberjack::log_status(string message){
    if (!logfile.is_open()){
        throw UnexpectedClosedFileException("Lumberjack::log_status");
    }
    logfile << get_time() << " status: " << message << endl;
}

void Lumberjack::log_debug(string message){
    if (!logfile.is_open()){
        throw UnexpectedClosedFileException("Lumberjack::log_debug");
    }
    logfile << get_time() << " status: " << message << endl;
}