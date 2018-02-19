
#ifndef SANA_ABSTRACT_HPP
#define SANA_ABSTRACT_HPP

#include <map>
#include "Method.hpp"
#include <random>
#include "../measures/localMeasures/LocalMeasure.hpp"
#include "../measures/Measure.hpp"
#include "../measures/MeasureCombination.hpp"
#include "../utils/randomSeed.hpp"
#include "../measures/ExternalWeightedEdgeConservation.hpp"

#include "../measures/SymmetricSubstructureScore.hpp"
#include "../measures/EdgeCorrectness.hpp"
#include "../measures/WeightedEdgeConservation.hpp"
#include "../measures/NodeCorrectness.hpp"
#include "../measures/SymmetricEdgeCoverage.hpp"
#include "../measures/localMeasures/Sequence.hpp"
#include "../utils/NormalDistribution.hpp"
#include "../utils/LinearRegression.hpp"
#include "../arguments/measureSelector.hpp"
#include "../measures/InducedConservedStructure.hpp"
#include "../measures/LargestCommonConnectedSubgraph.hpp"

#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <stdexcept>
#include <unordered_set>
#include <algorithm>
#include <random>
#include <queue>
#include <iomanip>
#include <set>
#include <cmath>
#include <cassert>
#include <signal.h>
#include <stdio.h>
#include <unistd.h>
#include <ctime>

class SANA_Abstract: public Method {
    public:
        SANA_Abstract(Graph* G1, Graph* G2,
                double TInitial, double TDecay, double t, bool usingIterations, bool addHillClimbing, MeasureCombination* MC, string& objectiveScore, string name);
        SANA_Abstract(const SANA_Abstract& orig);
        virtual ~SANA_Abstract();
        virtual Alignment run();
        virtual void describeParameters(ostream& stream);
        virtual string fileNameSuffix(const Alignment& A);
        
        virtual void enableRestartScheme(double minutesNewAlignments, uint iterationsPerStep,
                uint numCandidates, double minutesPerCandidate, double minutesFinalist);
        
        //set temperature schedule automatically
        virtual void setTemperatureScheduleAutomatically();
        virtual void setTInitialByLinearRegression();
        virtual void setTInitialByStatisticalTest();
        virtual void setTDecayAutomatically();
        
        //set temperature decay dynamically 
        virtual void setDynamicTDecay(); 
        
        double elapsedEstimate = 0;
        int order = 0;
        
        //returns the number of iterations until it stagnates when not using temperture
        virtual long long unsigned int hillClimbingIterations(long long unsigned int idleCountTarget);
        virtual Alignment hillClimbingAlignment(Alignment startAlignment, long long unsigned int idleCountTarget);
        virtual Alignment hillClimbingAlignment(long long unsigned int idleCountTarget);
        
        //returns an approximation of the the logarithm in base e of the size of the search space
        virtual double searchSpaceSizeLog();
        
    protected:    
        //Temperature Boundaries. Use these after the tinitial has been determined
        double lowerTBound = 0;
        double upperTBound = 0;

        //data structures for the networks
        uint n1;
        uint n2;
        double g1Edges; //stored as double because it appears in division
        double g2Edges; //stored as double because it appears in division
        vector<vector<bool> > G1AdjMatrix;
        vector<vector<bool> > G2AdjMatrix;
        vector<vector<ushort> > G1AdjLists;
        vector<vector<ushort> > G2AdjLists;

        virtual void initTau(void);
        vector<ushort> unLockedNodesG1;
        bool nodesHaveType = false;
        //random number generation
        mt19937 gen;
        uniform_int_distribution<> G1RandomNode;
        uniform_int_distribution<> G1RandomUnlockedNodeDist;
        uniform_int_distribution<> G2RandomUnassignedNode;
        uniform_real_distribution<> randomReal;
        uniform_int_distribution<> G1RandomUnlockedGeneDist;
        uniform_int_distribution<> G1RandomUnlockedmiRNADist;
        virtual ushort G1RandomUnlockedNode();
        virtual ushort G1RandomUnlockedNode(uint source1); // used in nodes-have-type because
        virtual ushort G1RandomUnlockedNode_Fast();
        virtual ushort G2RandomUnlockedNode(uint target1);
        virtual ushort G2RandomUnlockedNode_Fast();

        //temperature schedule
        double TInitial;
        double TDecay;
        double minutes = 0;
        bool usingIterations;
        uint maxIterations = 0;
        uint iterationsPerformed = 0;
        const double TInitialScaling = 1;
        const double TDecayScaling = 1;
        //to compute TDecay dynamically 
        //vector holds "ideal" temperature values at certain execution times 
        bool dynamic_tdecay;
        vector<double> tau; 
        double SANAtime;

        double T;
        virtual double temperatureFunction(double iter, double TInitial, double TDecay);
        virtual double acceptingProbability(double energyInc, double T);
        virtual double trueAcceptingProbability();
        //to compute TInitial automatically
        //returns a value of TInitial such that the temperature is random
        virtual double searchTInitialByStatisticalTest(), simpleSearchTInitial();
        virtual double scoreForTInitial(double TInitial);
        virtual bool isRandomTInitial(double TInitial, double highThresholdScore, double lowThresholdScore);
        virtual double scoreRandom();
        //to compute TDecay automatically
        //returns a value of lambda such that with this TInitial, temperature reaches
        //0 after a certain number of minutes
        virtual double searchTDecay(double TInitial, double minutes);
        virtual double searchTDecay(double TInitial, uint iterations);

        bool initializedIterPerSecond;
        double iterPerSecond;
        virtual double getIterPerSecond();
        virtual void initIterPerSecond();

        virtual vector<double> energyIncSample(double temp = 0.0);
        virtual double expectedNumAccEInc(double temp, const vector<double>& energyIncSample);

        //data structures for the solution space search
        double changeProbability;
        vector<bool> assignedNodesG2;
        vector<ushort> unassignedNodesG2;
        vector<ushort> A;
        //initializes all the necessary data structures for a new run
        virtual void initDataStructures(const Alignment& startA);

        bool addHillClimbing; //for post-run hill climbing

        //objective function
        MeasureCombination* MC;
        virtual double eval(const Alignment& A);
        virtual double scoreComparison(double newAligEdges, double newInducedEdges, double newLocalScoreSum, double newWecSum, double newNcSum, double& newCurrentScore, double newEwecSum);
        double ecWeight;
        double s3Weight;
        double wecWeight;
        double secWeight;
        double ncWeight;
        double localWeight;
        double ewecWeight;
        string score;


        //restart scheme
        bool restart;
        //parameters
        double minutesNewAlignments;
        uint iterationsPerStep;
        uint numCandidates;
        double minutesPerCandidate;
        double minutesFinalist;
        //data structures
        uint newAlignmentsCount;
        vector<Alignment> candidates;
        vector<double> candidatesScores;
        //functions
        virtual Alignment runRestartPhases();
        virtual uint getLowestIndex() const;
        virtual uint getHighestIndex() const;


        //to evaluate EC incrementally
        bool needAligEdges;
        int aligEdges;
        virtual int aligEdgesIncChangeOp(ushort source, ushort oldTarget, ushort newTarget);
        virtual int aligEdgesIncSwapOp(ushort source1, ushort source2, ushort target1, ushort target2);

        //to evaluate EC incrementally
        bool needSec;
        double secSum;

        //to evaluate S3 incrementally
        bool needInducedEdges;
        int inducedEdges;
        virtual int inducedEdgesIncChangeOp(ushort source, ushort oldTarget, ushort newTarget);

        //to evaluate nc incrementally
        bool needNC;
        int ncSum;
        vector<ushort> trueA;
        virtual int ncIncChangeOp(ushort source, ushort oldTarget, ushort newTarget);
        virtual int ncIncSwapOp(ushort source1, ushort source2, ushort target1, ushort target2);

        //to evaluate wec incrementally
        bool needWec;
        double wecSum;
        vector<vector<float> > wecSims;
        virtual double WECIncChangeOp(ushort source, ushort oldTarget, ushort newTarget);
        virtual double WECIncSwapOp(ushort source1, ushort source2, ushort target1, ushort target2);

        //to evaluate ewec incrementally
        bool needEwec;
        ExternalWeightedEdgeConservation* ewec;
        double ewecSum;
        virtual double EWECIncChangeOp(ushort source, ushort oldTarget, ushort newTarget, const Alignment& A);
        virtual double EWECIncSwapOp(ushort source1, ushort source2, ushort target1, ushort target2, const Alignment& A);

        //to evaluate local measures incrementally
        bool needLocal;
        double localScoreSum;
        map<string, double> localScoreSumMap;
        vector<vector<float> > sims;
        map<string, vector<vector<float> > > localSimMatrixMap;
        virtual double localScoreSumIncChangeOp(vector<vector<float> > const & sim, ushort const & source, ushort const & oldTarget, ushort const & newTarget);
        virtual double localScoreSumIncSwapOp(vector<vector<float> > const & sim, ushort const & source1, ushort const & source2, ushort const & target1, ushort const & target2);



        //other execution options
        bool constantTemp; //tempertare does not decrease as a function of iteration
        bool enableTrackProgress; //shows output periodically
        virtual void trackProgress(long long unsigned i);
        double avgEnergyInc;


        //algorithm
        virtual Alignment simpleRun(const Alignment& A, double maxExecutionSeconds,
                long long unsigned int& iter);
        virtual void evalAlignmentMethod();
        virtual vector<string> addFieldsToResults();

        virtual Alignment simpleRun(const Alignment& A, long long unsigned int maxExecutionIterations,
                long long unsigned int& iter);
        double currentScore;
        double energyInc;
        vector<double> sampledProbability;
        virtual void SANAIteration();
        virtual void performChange()=0;
        virtual void performSwap()=0;


        //others
        Timer timer;
        virtual void setInterruptSignal(); //allows the execution to be paused with Control+C

        // Used to support locking
        virtual Alignment getStartingAlignment()=0;
        virtual bool implementsLocking(){ return true; }

        virtual double pForTInitial(double TInitial);
        virtual double getPforTInitial(const Alignment& startA, double maxExecutionSeconds,
                long long unsigned int& iter);
        virtual double findTInitialByLinearRegression();
        virtual string getFolder();
        string haveFolder();
        virtual string mkdir(const std::string& file);
        //        tuple<int, double, int, double, double, double> regress(double start, double end, int amount);
            
};

#endif /* SANA_ABSTRACT_HPP */

