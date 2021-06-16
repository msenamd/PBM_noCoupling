#ifndef VEGREACTION_H
#define VEGREACTION_H

#include <string>
#include <iostream>
using namespace std;
#include <fstream>
#include <vector>
#include <tuple>
#include <algorithm>


class VegReaction{

public:
    VegReaction();
    ~VegReaction();
    //Reaction(string reactionName);


    bool set_knownReaction(std::string reactionName);

    double get_A();
    double get_Ta();
    double get_deltaH();

private:

    double A;
    double Ea;
    double deltaH;
};

#endif
