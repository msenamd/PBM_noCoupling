/*
 * reaction.h
 * A class of chemical reactions
 * All properties must be listed here
 * Different reaction mechanisms are to be implemented
 *
 *  Created on: October 7, 2019
 *      Author: mmahmed15
 */
#include "VegReaction.h"

VegReaction::VegReaction()
{

}

VegReaction::~VegReaction()
{

}


bool VegReaction::set_knownReaction(std::string reactionName)
{
    if(reactionName == "drying")  //  - Source: Shen et al. (2007) Fire Safety J. 42:210-217
    {
        A=5.13e+10;        //Pre-exponential factor [1/s]
        Ea=88e+3;          //Activation energy [J/mol]
        deltaH=-2440e+3;   //Heat of evaporation [J/kg]
        return true;
    }
    else if(reactionName == "pyrolysis") //- Source: Novozhilov et al. (1996) Fire Safety J. 27:69-84
    {
        A=5.25e+7;        //Pre-exponential factor [1/s]
        Ea=1.256e+5;      //Activation energy [J/mol]
        deltaH=-0.0;      //Heat of evaporation [J/kg
        return true;
    }
    else if(reactionName == "char oxidation") //- Source: Branca and Di Blasi (2004) Therm Sci 8:51-63
    {
        A=1.4e+11;         //Pre-exponential factor [1/s]
        Ea=180e+3;         //Activation energy [J/mol]
        deltaH=32e+6;      //Heat of evaporation [J/kg]
        return true;
    }
    else
    {
        cout << "Error: Invalid reaction" << endl;
        return false;   //unknown name
    }
}


double VegReaction::get_A()
{
    return A;
}

double VegReaction::get_Ta()
{
    return Ea/8.3145;
}

double VegReaction::get_deltaH()
{
    return deltaH;
}