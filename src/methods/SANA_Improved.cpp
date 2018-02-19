
#include "SANA_Improved.hpp"

using namespace std;

void SANA_Improved::initBucketsDataStructurs(){
    vector<uint> degreeVector = getDegreeDistribution(G2);
    map<int, vector<ushort>> degreeToNode;    
    for (uint i = 0; i < degreeVector.size(); i++) {
        degreeToNode[degreeVector[i]].push_back(i);
    } 
    if(NULL == buckets)
        delete(buckets);
    if(NULL == bucket_sizes)
        delete(bucket_sizes);
    if(NULL == node_to_bucket)
        delete(node_to_bucket);
    buckets = new vector<vector<ushort>>();
    bucket_sizes = new vector<vector<int>>();    
    node_to_bucket = new vector<vector<int>>(n2, vector<int>(2));
    map<int,vector<ushort> >::iterator deg;
    int b = 0;
    for(deg = degreeToNode.begin(); deg != degreeToNode.end(); deg++)
    {
        for(uint node_index = 0; node_index < deg->second.size(); node_index++)
            node_to_bucket->at(deg->second[node_index]) = {b, node_index};
        buckets->push_back(deg->second);
        bucket_sizes->push_back({0, (int)deg->second.size()});
        b++;
    }
}

void SANA_Improved::initBucketAccSize(){
    vector<vector<int>> t_bucket_sizes_acc;
    int free_sum = 0;
    int non_free_sum = 0;
    for(uint j = 0; j < bucket_sizes->size(); j++)
    {
        t_bucket_sizes_acc.push_back({non_free_sum, free_sum});
        non_free_sum += bucket_sizes->at(j)[0];
        free_sum += bucket_sizes->at(j)[1];
    }
    bucket_sizes_acc = new vector<vector<int>>(t_bucket_sizes_acc); 
}

void SANA_Improved::bucket_extend_fn(map<int,vector<int> >& degreeToBucket){
   vector<uint> g2DegreeVector = getDegreeDistribution(G2);
    for(map<int, vector<int>>::iterator i = degreeToBucket.begin(); i != degreeToBucket.end(); i++)
    {
        int min_b = i->second[0];
        for(int j = min_b - 1; j >= max(0, (min_b - sqrt(min_b))); j--)
            i->second.insert(i->second.begin(), j);
        int max_b = i->second[i->second.size() - 1];
        int j = max_b + 1;
        int max_b_deg = g2DegreeVector[buckets->at(max_b)[0]];
	max_b++;
        while(i->first > max_b_deg && bucket_sizes->size() > max_b){            
            max_b_deg = g2DegreeVector[buckets->at(max_b)[0]];
	    max_b++;
        }
        for(; j <= min(bucket_sizes->size() - 1, (max_b + sqrt(max_b))); j++)
            i->second.push_back(j);
    }
}

void SANA_Improved::initAlignment(Graph* G1, Graph* G2, vector<ushort>& alignment, vector<ushort>& n_alignment){
    if(NULL == candidates_for_node)
        delete(candidates_for_node);
    candidates_for_node = new vector<vector<int>>(n1);
    map<int, vector<int>> degreeToBucket;
    vector<uint> g1SortedNodes = getNodesSortedByDegree(G1);
    vector<uint> g2SortedNodes = getNodesSortedByDegree(G2);       
    vector<uint> g1DegreeVector = getDegreeDistribution(G1);
    int curr_degree = -1;
    int curr_bucket = -1;
    for(uint i = 0; i < g1SortedNodes.size(); i++)
    {        
        uint v1 = g1SortedNodes[i];
        uint v2 = g2SortedNodes[i];
        alignment[v1] = v2;
        n_alignment[v2] = v1;
        int b = node_to_bucket->at(v2)[0];
        moveToNonFree(v2);
        int d = g1DegreeVector[g1SortedNodes[i]];
        if(curr_bucket != b || curr_degree != d)
            degreeToBucket[g1DegreeVector[v1]].push_back(b);
        curr_bucket = b;
        curr_degree = d;        
    }    
    

    bucket_extend_fn(degreeToBucket);
    
    map<int, vector<int>> degreeToNode;
    for(int i = 0; i < n1; i++)
    {
        int deg = g1DegreeVector[i];
        degreeToNode[deg].push_back(i);
    }
    
    
    ofstream deg_to_b;
    deg_to_b.open("DTB_"+ G1->getName()+"_"+G2->getName()+".txt");
    for(map<int, vector<int>>::iterator i = degreeToBucket.begin(); i != degreeToBucket.end(); i++)
    {
       int deg = i->first;
       deg_to_b << "Deg: " << deg<<" Deg size: "<< float(degreeToNode[deg].size()) / n1 * 100 << "% Buckets Size: ";
       int count = 0;
       for(int j = 0; j < i->second.size(); j++)
       {
           count += buckets->at(i->second[j]).size();
       }
       deg_to_b << float(count) / n2 * 100 <<"%"<<endl;
    }
    deg_to_b.flush();
    deg_to_b.close();
    
    
    for(uint node = 0; node < n1; node++)
    {
        candidates_for_node->at(node) = degreeToBucket[g1DegreeVector[node]];
    }
}

SANA_Improved::SANA_Improved(Graph* G1, Graph* G2,
        double TInitial, double TDecay, double t, bool usingIterations, bool addHillClimbing, MeasureCombination* MC, string& objectiveScore):
SANA_Abstract(G1, G2, TInitial, TDecay, t, usingIterations, addHillClimbing, MC, objectiveScore, "SANA_Improved")
{
    double startTime = timer.elapsed();
    
    //init alignment
    
    initBucketsDataStructurs();
    vector<ushort> alignment(n1);   
    vector<ushort> n_alignment(n2);     
    initAlignment(G1, G2, alignment, n_alignment);
    curr_alignment = alignment;
    initBucketAccSize();
    
    
    vector<string> g1NodeNames = G1->getNodeNames();
    map<string,ushort> g1NodeToIndexMap = G1->getNodeNameToIndexMap();
    
    vector<string> g2NodeNames = G2->getNodeNames();
    map<string,ushort> g2NodeToIndexMap = G2->getNodeNameToIndexMap();
    
    vector<uint> g2DegreeVector = getDegreeDistribution(G2);
    
    igraph_vector_t g1Cores;
    igraph_t graph1;          
    initIgraphFromGraph(&graph1, G1);    
    igraph_vector_init(&g1Cores, 0);
    igraph_coreness(&graph1, &g1Cores, IGRAPH_ALL);
    
    igraph_vector_t g2Cores;
    igraph_t graph2;          
    initIgraphFromGraph(&graph2, G2);    
    igraph_vector_init(&g2Cores, 0);
    igraph_coreness(&graph2, &g2Cores, IGRAPH_ALL);   
    vector<uint> g1DegreeVector = getDegreeDistribution(G1);
    
       
    double endTime = timer.elapsed();
    cout << "***** preprossing time: " << endTime - startTime <<endl; 
}

template<typename T> void safe_delete(T*& a) {
  delete a;
  a = NULL;
}

SANA_Improved::~SANA_Improved() {      
    safe_delete(candidates_for_node);    
    safe_delete(buckets);    
    safe_delete(bucket_sizes);    
    safe_delete(bucket_sizes_acc);    
    safe_delete(node_to_bucket);    
}


void SANA_Improved::initDataStructures(const Alignment& startA)
{
    NA =  vector<int>(n2);
    Alignment startA_temp = startA;
    SANA_Abstract::initDataStructures(startA_temp);
    for(ushort node = 0; node < n1; node++)
    {
        int u1 = A[node];
        NA[u1] = node;
    }    
}

vector<uint> SANA_Improved::getDegreeDistribution(Graph* G)
{
    igraph_t graph;      
    vector<uint> resultVector;
    initIgraphFromGraph(&graph, G);
    igraph_vector_t degreeVector;  
    igraph_vector_init(&degreeVector, 0);
    igraph_degree(&graph, &degreeVector, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);  
    for (uint i=0; i<igraph_vector_size(&degreeVector); i++) {
        resultVector.push_back(VECTOR(degreeVector)[i]);
    } 
    return resultVector;
}

vector<uint> SANA_Improved::getNodesSortedByDegree(Graph* G) 
{    
    vector<uint> resultVector;
    vector<uint> degreeVector = getDegreeDistribution(G);
    map<int, vector<uint>> degreeToNode;
    
    for (uint i = 0; i < degreeVector.size(); i++) {
        degreeToNode[degreeVector[i]].push_back(i);
    } 
    map<int,vector<uint> >::iterator i;
    for(i = degreeToNode.begin(); i != degreeToNode.end(); i++)
    {
        for(uint j = 0; j < i->second.size(); j++)
        {
            resultVector.push_back(i->second[j]);
        }
    }
    return resultVector;
}

void SANA_Improved::initIgraphFromGraph(igraph_t * graph, Graph* G)
{
    igraph_vector_t edges_vector;
    int vertices_count = G->getNumNodes();
    initIgraphEdgeVectorFromGraph(&edges_vector, G);
    
    igraph_empty(graph, vertices_count, IGRAPH_UNDIRECTED);
    
    igraph_add_edges(graph, &edges_vector, 0);
    
}

void SANA_Improved::initIgraphEdgeVectorFromGraph(igraph_vector_t *edges_vector, Graph* G)
{
    int edgesCount = G->getNumEdges();
    vector<vector<ushort> >edges (edgesCount);
    vector<vector<ushort> > &edgeList = edges;
    G->getEdgeList(edgeList);
    
    igraph_vector_init(edges_vector, edgesCount*2);
    int i = 0;
    //  convert edgeList to edges vector for igraph
    for (vector<vector<ushort>>::iterator it = edgeList.begin() ; it != edgeList.end(); ++it)
    {        
        VECTOR(*edges_vector)[i++] = (*it)[0];
        VECTOR(*edges_vector)[i++] = (*it)[1];
    }
}

Alignment SANA_Improved::getStartingAlignment(){  
    initBucketsDataStructurs();
    vector<ushort> alignment(n1);   
    vector<ushort> n_alignment(n2);     
    initAlignment(G1, G2, alignment, n_alignment);
    curr_alignment = alignment;
    initBucketAccSize();
        
    return alignment;
}

void SANA_Improved::moveToNonFree(ushort node){
    int b_id = node_to_bucket->at(node)[0];
    int node_index_in_b = node_to_bucket->at(node)[1];
    
    int& nonFreeSize = bucket_sizes->at(b_id)[0];
    int& freeSize = bucket_sizes->at(b_id)[1];        
    
    vector<ushort>& bucket_nodes = buckets->at(b_id);
    ushort temp = bucket_nodes[nonFreeSize];
    bucket_nodes[node_index_in_b] = temp;
    bucket_nodes[nonFreeSize] = node;
    node_to_bucket->at(node)[1] = nonFreeSize;
    node_to_bucket->at(temp)[1] = node_index_in_b;
    nonFreeSize++;
    freeSize--;
}

void SANA_Improved::moveToFree(ushort node){
    int b_id = node_to_bucket->at(node)[0];
    int node_index_in_b = node_to_bucket->at(node)[1];
    
    int& nonFreeSize = bucket_sizes->at(b_id)[0];
    int& freeSize = bucket_sizes->at(b_id)[1];  
    
    if(nonFreeSize != 0){
        vector<ushort>& bucket_nodes = buckets->at(b_id);
        int lastNonFree = nonFreeSize - 1;
        ushort temp = bucket_nodes[lastNonFree];
        bucket_nodes[node_index_in_b] = temp;
        bucket_nodes[lastNonFree] = node;
        node_to_bucket->at(node)[1] = lastNonFree;
        node_to_bucket->at(temp)[1] = node_index_in_b;
        nonFreeSize--;
        freeSize++;
    }    
}

void SANA_Improved::printBucketAcc(){
    for(int b = 0; b < bucket_sizes_acc->size(); b++)
    {
        int acc_non = bucket_sizes_acc->at(b)[0];
        int acc_free = bucket_sizes_acc->at(b)[1];
        cout << "Bucket: "<< b << " acc non: " << acc_non;
        cout << " acc free: " << acc_free << endl;
    }
}

void SANA_Improved::validate_alignment(ushort v1, ushort u1){
    if(NA[u1] != v1){
        cout << "Error" <<endl;
        cout << "NA[u1]: " << NA[u1] << " v1: " << v1 << endl;
        exit(-1);
    }
}

void SANA_Improved::validate_A(){
    int count = 0;
    for(int i = 0; i< A.size(); i++)
    {
        ushort j = A[i];
        if(NA[j] != i)
            count++;
    }
    if(count > 0)
        cout << "Not good count: " << count << endl;
}

void SANA_Improved::performChange() {
    static int count = 0;
    count++;
    ushort v1 = G1RandomUnlockedNode();
    ushort u1 = A[v1];
    
       
    vector<int> candidates1 = candidates_for_node->at(v1);
    int start_b_v1 = candidates1[0];        
    int end_b_v1 = candidates1[candidates1.size() - 1];  
    
    vector<int> newTarget = getRandomBucketAndIndexOfNode(start_b_v1, end_b_v1, 1); 
    if(newTarget[0] != -1){
    
        int u2_b_id = newTarget[0];
        int u2_index = newTarget[1];
        
        ushort u2 =  buckets->at(u2_b_id)[u2_index];  
        
        validate_alignment(v1, u1);        
    
        int newAligEdges = -1; //dummy initialization to shut compiler warnings
        if (needAligEdges or needSec) {
            newAligEdges = aligEdges + aligEdgesIncChangeOp(v1, u1, u2);
        }

        int newInducedEdges = -1; //dummy initialization to shut compiler warnings
        if (needInducedEdges) {
            newInducedEdges = inducedEdges + inducedEdgesIncChangeOp(v1, u1, u2);
        }

        double newLocalScoreSum = -1; //dummy initialization to shut compiler warnings
        map<string, double> newLocalScoreSumMap(localScoreSumMap);
        if (needLocal) {
            newLocalScoreSum = localScoreSum + localScoreSumIncChangeOp(sims, v1, u1, u2);
            for(auto it = newLocalScoreSumMap.begin(); it != newLocalScoreSumMap.end(); ++it)
                it->second += localScoreSumIncChangeOp(localSimMatrixMap[it->first], v1, u1, u2);
        }

        double newWecSum = -1; //dummy initialization to shut compiler warning
        if (needWec) {
            newWecSum = wecSum + WECIncChangeOp(v1, u1, u2);
        }

        double newEwecSum = -1;
        if (needEwec) {
            newEwecSum = ewecSum + EWECIncChangeOp(v1, u1, u2, A);
        }

        double newNcSum = -1;
        if (needNC) {
            newNcSum = ncSum + ncIncChangeOp(v1, u1, u2);
        }	

        double newCurrentScore = 0;
        bool makeChange = scoreComparison(newAligEdges, newInducedEdges, newLocalScoreSum, newWecSum, newNcSum, newCurrentScore, newEwecSum);

        if (makeChange) {
            A[v1] = u2;  
            NA[u2] = v1;
            NA[u1] = -1;
                        
            moveToNonFree(u2); 
            moveToFree(u1);     
          
            int b1 = node_to_bucket->at(u1)[0];
            int b2 = node_to_bucket->at(u2)[0];
            
            if (b1 < b2){
                for(int b = b1 + 1; b <= b2; b++){
                    bucket_sizes_acc->at(b)[0]--;
                    bucket_sizes_acc->at(b)[1]++;
                }
            }
            else{
                for(int b = b2 + 1; b <= b1; b++){
                    bucket_sizes_acc->at(b)[0]++;
                    bucket_sizes_acc->at(b)[1]--;
                }
            }

            assignedNodesG2[u1] = false;
            assignedNodesG2[u2] = true;
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
}

void SANA_Improved::performSwap() { 
    
    static int count = 0;
    count++;
    ushort v1 = G1RandomUnlockedNode();   
    ushort u1 =  A[v1];   
// 
    vector<int> candidates1 = candidates_for_node->at(v1);
    int start_b_v1 = candidates1[0];        
    int end_b_v1 = candidates1[candidates1.size() - 1];  
    vector<int> newTarget = getRandomBucketAndIndexOfNode(start_b_v1, end_b_v1, 0);          
    
    int bNew = newTarget[0];
    int indexNew = newTarget[1];  
    ushort u2 =  buckets->at(bNew)[indexNew];    
    ushort v2 = NA[u2];

    vector<int> candidates2 = candidates_for_node->at(v2);
    int u1_bucket = node_to_bucket->at(u1)[0];
    int start_candidate_v2 = candidates2[candidates2.size() - 1];    
    int end_candidate_v2 = candidates2[0];    
    if(u1_bucket <= start_candidate_v2 || u1_bucket >= end_candidate_v2)
    {  
        int newAligEdges = -1; //dummy initialization to shut compiler warnings
        if (needAligEdges or needSec) {
            newAligEdges = aligEdges + aligEdgesIncSwapOp(v1, v2, u1, u2);
        }
        
        double newLocalScoreSum = -1; //dummy initialization to shut compiler warnings
        map<string, double> newLocalScoreSumMap(localScoreSumMap);
        if (needLocal) {
            newLocalScoreSum = localScoreSum + localScoreSumIncSwapOp(sims, v1, v2, u1, u2);
            for(auto it = newLocalScoreSumMap.begin(); it != newLocalScoreSumMap.end(); ++it)
                it->second += localScoreSumIncSwapOp(localSimMatrixMap[it->first], v1, v2, u1, u2);
        }
        
        double newWecSum = -1; //dummy initialization to shut compiler warning
        if (needWec) {
            newWecSum = wecSum + WECIncSwapOp(v1, v2, u1, u2);
        }
        
        double newEwecSum = -1;
        if (needEwec) {
            newEwecSum = ewecSum + EWECIncSwapOp(v1, v2, u1, u2, A);
        }
        
        
        double newNcSum = -1;
        if(needNC) {
            newNcSum = ncSum + ncIncSwapOp(v1, v2, u1, u2);
        }
        double newCurrentScore = 0;
        bool makeChange = scoreComparison(newAligEdges, inducedEdges, newLocalScoreSum, newWecSum, newNcSum, newCurrentScore, newEwecSum);
        
        if (makeChange) {
            A[v1] = u2;
            A[v2] = u1;
            
            NA[u2] = v1;
            NA[u1] = v2;
            
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
}


vector<int> SANA_Improved::getRandomBucketAndIndexOfNode(int start_bucket, int end_bucket, int kind)
{
    int end_size_acc = bucket_sizes_acc->at(end_bucket)[kind];
    int start_size_acc = bucket_sizes_acc->at(start_bucket)[kind];
    int end_size = bucket_sizes->at(end_bucket)[kind];
    int total = end_size_acc - start_size_acc + end_size;    
    if(total == 0)
        return {-1, -1};
    int rand_index = start_size_acc + rand() % total;   
    int b = searchForBucketOfIndex(start_bucket, end_bucket, kind, rand_index);   
    int node_index;
    int posible_index = rand_index - bucket_sizes_acc->at(b)[kind];
    if(kind == 0) 
    {
        node_index = posible_index;       
    }
    else{
        int first_free_index = bucket_sizes->at(b)[0];
        node_index = posible_index + first_free_index;
    }
    vector<int> result_bucket_node_index = {b, node_index};    
    return result_bucket_node_index;
}

int SANA_Improved::searchForBucketOfIndex(int start_bucket, int end_bucket, int kind, int rand_index){
    int mid_bucket = (start_bucket + end_bucket) / 2;
    int mid_bucket_size = bucket_sizes->at(mid_bucket)[kind];
    int posible_index = rand_index - bucket_sizes_acc->at(mid_bucket)[kind];    
    if(posible_index < 0){
        searchForBucketOfIndex(start_bucket, mid_bucket - 1, kind, rand_index);
    }
    else if(posible_index < mid_bucket_size)
    {        
        
        return mid_bucket;
    }
    else{
        searchForBucketOfIndex(mid_bucket + 1, end_bucket, kind, rand_index);
    }
}

inline bool isFileExists (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

void SANA_Improved::evalAlignmentMethod(){
    SANA_Abstract::evalAlignmentMethod();
}

vector<string> SANA_Improved::addFieldsToResults()
{ 
    map<ushort,string> G1Index2Name = G1->getIndexToNodeNameMap();
    map<string,ushort> G2Name2Index = G2->getNodeNameToIndexMap();

    uint n = G1->getNumNodes();
    vector<ushort> correctAlignment(n);
    for (uint i = 0; i < n; i++) {
        correctAlignment[i] = G2Name2Index[G1Index2Name[i]];
    }
    float OP;
    float in_bucket = 0;
    for(int node = 0; node < n1; node++)
    {
        int alignNode = correctAlignment[node];
        vector<int> candidates = candidates_for_node->at(node);
        int minB = candidates[0];
        int maxB = candidates[candidates.size() - 1];
        int alignNodeB = node_to_bucket->at(alignNode)[0];
        if(alignNodeB >= minB && alignNodeB <= maxB)
            in_bucket++;
    }
    OP = (in_bucket / n1) * 100;
    string op = to_string(OP) + "%";
    vector<string> fields;
    fields.push_back(op);
    return fields;
}


