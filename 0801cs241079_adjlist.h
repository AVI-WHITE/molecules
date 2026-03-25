#include<stdio.h>
#include<stdlib.h>
#define SPARE 5

// Node of the Linked list
struct LinkNode
{
    int vertexIndex;
    int edgeWeight;
    struct LinkNode *next;
};
typedef struct LinkNode SllNode;

typedef struct AdjacencyList{
    // array of SllNodes
    SllNode* arr;
    // number of nodes in the arr
    int nodeCount;
    // total capacity of array
    int capacity;
    // type of graph : 0 for directed and 1 for undirected
    int graphType;
} AdjacencyList;    

// function declarations of linked list functions
SllNode* addNodeBegin(SllNode* head, int vi, int wt);
SllNode* addSllNodeAfterGivenNode(SllNode* head, int vi, int wt);
SllNode* createSllNode(int vi, int wt);
  
// type = 1 for undirected graphs and 0 for directed graphs
AdjacencyList* createNewAdjacencyList(int numVertices, int type){
    AdjacencyList* graph = (AdjacencyList *)malloc(sizeof(AdjacencyList));

    graph->arr = (SllNode *)malloc(sizeof(SllNode) * (numVertices + SPARE));
    graph->capacity = numVertices + SPARE;
    graph->nodeCount = numVertices;
    graph->graphType = type;

    // setting the next of all the nodes to NULL
    for(int i=0; i < graph->capacity; i++){
        graph->arr[i].vertexIndex = i;
        graph->arr[i].next = NULL;
    }

    return graph;
}    

void addEdgeAL(AdjacencyList* graph, int v1, int v2, int weight){
    if(v1 >= graph->nodeCount || v2 >= graph->nodeCount || v1 == v2){
        return;
    }
    graph->arr[v1].next = addNodeBegin(graph->arr[v1].next, v2, weight);
    if(graph->graphType){
        graph->arr[v2].next = addNodeBegin(graph->arr[v2].next, v1, weight);
    }
}

int addVertexAL(AdjacencyList* graph){
    if(graph->nodeCount == graph->capacity){
        graph->capacity += SPARE;
        SllNode* newArr = (SllNode *)realloc(graph->arr, graph->capacity * sizeof(SllNode));
        graph->arr = newArr;
    }

    // adding new vertex
    graph->arr[graph->nodeCount].vertexIndex = graph->nodeCount;
    graph->arr[graph->nodeCount].next = NULL;
    return graph->nodeCount++;
}

void printAdjacencyList(AdjacencyList* graph){
    SllNode* temp;
    for(int i = 0; i < graph->nodeCount; i++){
        printf("\n%3d-> ", i);
        temp = graph->arr[i].next;
        while(temp != NULL){
            printf("(%3d, %3d) ", temp->vertexIndex, temp->edgeWeight);
            temp = temp->next;
        }
    }
}

// functions for linked list
SllNode* createSllNode(int vi, int wt)
{
    SllNode* newNode = (SllNode*)malloc(sizeof(SllNode));
    newNode->vertexIndex = vi;
    newNode->edgeWeight = wt;
    newNode->next = NULL;
    return newNode;
}

SllNode* addSllNodeAfterGivenNode(SllNode* head, int vi, int wt)
{
    if(head == NULL)
        return NULL;
    
    SllNode* newNode = createSllNode(vi, wt);
    newNode->next = head->next;
    head->next = newNode;
    return head;
}

SllNode* addNodeBegin(SllNode* head, int vi, int wt)
{
    SllNode* newNode = createSllNode(vi, wt);
    newNode->next = head;
    return newNode;
}