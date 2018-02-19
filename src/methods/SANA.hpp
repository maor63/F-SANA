#ifndef SANA_HPP
#define SANA_HPP
#include "SANA_Abstract.hpp"

class SANA: public SANA_Abstract {

public:
    SANA(Graph* G1, Graph* G2,
        double TInitial, double TDecay, double t, bool usingIterations, bool addHillClimbing, MeasureCombination* MC, string& objectiveScore);
    ~SANA();
    
protected: 
    void performChange();
    void performSwap();
    string haveFolder();
    void evalAlignmentMethod();
    virtual vector<string> addFieldsToResults();
    
    // Used to support locking
    virtual Alignment getStartingAlignment();
    
};

#endif
