

#include "SANA.hpp"

using namespace std;

SANA::SANA(Graph* G1, Graph* G2,
        double TInitial, double TDecay, double t, bool usingIterations, bool addHillClimbing, MeasureCombination* MC, string& objectiveScore):
SANA_Abstract(G1, G2, TInitial, TDecay, t, usingIterations, addHillClimbing, MC, objectiveScore, "SANA")
{
    
}

SANA::~SANA() {
    
}


Alignment SANA::getStartingAlignment(){
    if(G1->hasNodeTypes()){
        Alignment randomAlig = Alignment::randomAlignmentWithNodeType(G1,G2);
        randomAlig.reIndexBefore_Iterations(G1->getNodeTypes_ReIndexMap());
        return randomAlig;
    }
    else if(lockFileName != ""){
        Alignment randomAlig = Alignment::randomAlignmentWithLocking(G1,G2);
        randomAlig.reIndexBefore_Iterations(G1->getLocking_ReIndexMap());
        return randomAlig;
    }
    else{
        return Alignment::random(n1, n2);
    }
}

void SANA::performChange() {
    ushort source = G1RandomUnlockedNode();
    ushort oldTarget = A[source];
    
    uint newTargetIndex =  G2RandomUnlockedNode(oldTarget);
    ushort newTarget = unassignedNodesG2[newTargetIndex];
    
    //	assert(!G1->isLocked(source));
    //	assert(!G2->isLocked(newTarget));
    
    //	bool G1Gene = source < G1->unlockedGeneCount;
    //	bool G2Gene =  G2->nodeTypes[newTarget] == "gene";
    //	assert((G1Gene && G2Gene) || (!G1Gene && !G2Gene));
    
    
    int newAligEdges = -1; //dummy initialization to shut compiler warnings
    if (needAligEdges or needSec) {
        newAligEdges = aligEdges + aligEdgesIncChangeOp(source, oldTarget, newTarget);
    }
    
    int newInducedEdges = -1; //dummy initialization to shut compiler warnings
    if (needInducedEdges) {
        newInducedEdges = inducedEdges + inducedEdgesIncChangeOp(source, oldTarget, newTarget);
    }
    
    double newLocalScoreSum = -1; //dummy initialization to shut compiler warnings
    map<string, double> newLocalScoreSumMap(localScoreSumMap);
    if (needLocal) {
        newLocalScoreSum = localScoreSum + localScoreSumIncChangeOp(sims, source, oldTarget, newTarget);
        for(auto it = newLocalScoreSumMap.begin(); it != newLocalScoreSumMap.end(); ++it)
            it->second += localScoreSumIncChangeOp(localSimMatrixMap[it->first], source, oldTarget, newTarget);
    }
    
    double newWecSum = -1; //dummy initialization to shut compiler warning
    if (needWec) {
        newWecSum = wecSum + WECIncChangeOp(source, oldTarget, newTarget);
    }
    
    double newEwecSum = -1;
    if (needEwec) {
        newEwecSum = ewecSum + EWECIncChangeOp(source, oldTarget, newTarget, A);
    }
    
    double newNcSum = -1;
    if (needNC) {
        newNcSum = ncSum + ncIncChangeOp(source, oldTarget, newTarget);
    }	
    
    double newCurrentScore = 0;
    bool makeChange = scoreComparison(newAligEdges, newInducedEdges, newLocalScoreSum, newWecSum, newNcSum, newCurrentScore, newEwecSum);
    
    if (makeChange) {
        A[source] = newTarget;
        unassignedNodesG2[newTargetIndex] = oldTarget;
        assignedNodesG2[oldTarget] = false;
        assignedNodesG2[newTarget] = true;
        aligEdges = newAligEdges;
        inducedEdges = newInducedEdges;
        localScoreSum = newLocalScoreSum;
        for(auto const & newLocalScoreSumEntry : newLocalScoreSumMap)
            localScoreSumMap[newLocalScoreSumEntry.first] = newLocalScoreSumEntry.second;
        wecSum = newWecSum;
        ewecSum = newEwecSum;
        ncSum = newNcSum;
        currentScore = newCurrentScore;
    }
}

void SANA::performSwap() {
    ushort source1 =  G1RandomUnlockedNode();
    ushort source2 =  G1RandomUnlockedNode(source1);
    ushort target1 = A[source1], target2 = A[source2];
    
    //	if(!(source1 >= 0 and source1 < G1->getNumNodes()) || !(source2 >= 0 and source2 < G1->getNumNodes())){
    //	    cerr << source1 << "   " << source2 << endl;
    //	    cerr << G1->getNumNodes() << endl;
    //	}
    //	assert(source1 >= 0 and source1 < G1->getNumNodes());
    //	assert(source2 >= 0 and source2 < G1->getNumNodes());
    //
    //	bool s1Gene = source1 < G1->unlockedGeneCount;
    //    bool s2Gene = source2 < G1->unlockedGeneCount;
    //    if(not((s1Gene && s2Gene) || (!s1Gene && !s2Gene))){
    //        cerr << source1 << " " << source2 << endl;
    //        cerr << G1->unlockedGeneCount << "   " << G1->getNumNodes() << endl;
    //    }
    //    assert((s1Gene && s2Gene) || (!s1Gene && !s2Gene));
    
    
    int newAligEdges = -1; //dummy initialization to shut compiler warnings
    if (needAligEdges or needSec) {
        newAligEdges = aligEdges + aligEdgesIncSwapOp(source1, source2, target1, target2);
    }
    
    double newLocalScoreSum = -1; //dummy initialization to shut compiler warnings
    map<string, double> newLocalScoreSumMap(localScoreSumMap);
    if (needLocal) {
        newLocalScoreSum = localScoreSum + localScoreSumIncSwapOp(sims, source1, source2, target1, target2);
        for(auto it = newLocalScoreSumMap.begin(); it != newLocalScoreSumMap.end(); ++it)
            it->second += localScoreSumIncSwapOp(localSimMatrixMap[it->first], source1, source2, target1, target2);
    }
    
    double newWecSum = -1; //dummy initialization to shut compiler warning
    if (needWec) {
        newWecSum = wecSum + WECIncSwapOp(source1, source2, target1, target2);
    }
    
    double newEwecSum = -1;
    if (needEwec) {
        newEwecSum = ewecSum + EWECIncSwapOp(source1, source2, target1, target2, A);
    }
    
    
    double newNcSum = -1;
    if(needNC) {
        newNcSum = ncSum + ncIncSwapOp(source1, source2, target1, target2);
    }
    double newCurrentScore = 0;
    bool makeChange = scoreComparison(newAligEdges, inducedEdges, newLocalScoreSum, newWecSum, newNcSum, newCurrentScore, newEwecSum);
    
    if (makeChange) {
        A[source1] = target2;
        A[source2] = target1;
        aligEdges = newAligEdges;
        localScoreSum = newLocalScoreSum;
        for(auto const & newLocalScoreSumEntry : newLocalScoreSumMap)
            localScoreSumMap[newLocalScoreSumEntry.first] = newLocalScoreSumEntry.second;
        wecSum = newWecSum;
        ewecSum = newEwecSum;
        ncSum = newNcSum;
        currentScore = newCurrentScore;
    }
}

void SANA::evalAlignmentMethod(){
    SANA_Abstract::evalAlignmentMethod();
}

vector<string> SANA::addFieldsToResults()
{
    vector<string> fields;
    string op = "100%";
    fields.push_back(op);
    return fields;
}

