

#ifndef SANA_IMPROVED_HPP
#define SANA_IMPROVED_HPP

#include "SANA_Abstract.hpp"
#include <igraph/igraph.h>

class SANA_Improved: public SANA_Abstract {
enum Flag { UNASSIGNED = 0, ASSIGNED = 1};

public:
    SANA_Improved(Graph* G1, Graph* G2,
        double TInitial, double TDecay, double t, bool usingIterations, bool addHillClimbing, MeasureCombination* MC, string& objectiveScore);
    void initAlignment(Graph* G1, Graph* G2, vector<ushort>& alignment, vector<ushort>& n_alignment);
    virtual void bucket_extend_fn(map<int,vector<int> >& degreeToBucket);


    void initBucketAccSize();

    void initBucketsDataStructurs();
    
    void initBucketsDataStructursOneBucket();

    ~SANA_Improved();
    
    

protected: 
    void performChange();
    void validate_A();

    void validate_alignment(ushort v1, ushort u1);

    void printBucketAcc();

    void moveToFree(ushort node);

    void moveToNonFree(ushort node);

    void performSwap();
    
    void evalAlignmentMethod();
    virtual vector<string> addFieldsToResults();
    
    string haveFolder();
    
    // Used to support locking
    virtual Alignment getStartingAlignment();
    vector<uint> getNodesSortedByDegree(Graph* G);
    void initIgraphFromGraph(igraph_t * graph, Graph* G);
    void initIgraphEdgeVectorFromGraph(igraph_vector_t *edges_vector, Graph* G);
    virtual void initDataStructures(const Alignment& startA);   
    vector<uint> getDegreeDistribution(Graph* G);
    vector<int> getRandomBucketAndIndexOfNode(int start_bucket, int end_bucket, int kind);
    int searchForBucketOfIndex(int start_bucket, int end_bucket, int kind, int rand_index);
        
    vector<ushort> curr_alignment;
    vector<vector<int>> *candidates_for_node; //vector
    vector<vector<ushort>> *buckets; //vector
    vector<vector<int>> *bucket_sizes; //vector
    vector<vector<int>> *bucket_sizes_acc; //vector
    vector<vector<int>> *node_to_bucket; //vector
    vector<int> NA;   
};

#endif /* SANA_IMPROVED_HPP */

