#pragma once

#include<stdio.h>
#include<stdlib.h>
#define SPARE 5

// Node of the Linked list
struct LinkNode
{
    int vertexIndex;
    int bondID;
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
    // total edges in graph
    int bondCount;
    // total capacity of bonds
    int bondCapacity;

} AdjacencyList;    

// function declarations of linked list functions
SllNode* addNodeBegin(SllNode* head, int vi, int bi);
SllNode* createSllNode(int vi, int bi);
void deleteLL(SllNode* head);
  

AdjacencyList* createNewAdjacencyList(int numVertices){
    AdjacencyList* graph = (AdjacencyList *)malloc(sizeof(AdjacencyList));

    graph->arr = (SllNode *)malloc(sizeof(SllNode) * (numVertices + SPARE));
    graph->capacity = numVertices + SPARE;
    graph->nodeCount = numVertices;
    graph->bondCount = 0;
    graph->bondCapacity = SPARE;

    // setting the next of all the nodes to NULL
    for(int i=0; i < graph->capacity; i++){
        graph->arr[i].vertexIndex = i;
        graph->arr[i].next = NULL;
    }

    return graph;
}    

int addEdgeAL(AdjacencyList* graph, int v1, int v2){
    if(v1 >= graph->nodeCount || v2 >= graph->nodeCount || v1 == v2){
        return -1;
    }
    graph->arr[v1].next = addNodeBegin(graph->arr[v1].next, v2, graph->bondCount);

    graph->arr[v2].next = addNodeBegin(graph->arr[v2].next, v1, graph->bondCount);

    return graph->bondCount++;
}

int addVertexAL(AdjacencyList* graph){
    if(graph->nodeCount == graph->capacity){
        graph->capacity += SPARE;
        SllNode* newArr = (SllNode *)realloc(graph->arr, graph->capacity * sizeof(SllNode));

        if(newArr == NULL){
            printf("\nMemory Allocation Failed!");
            exit(1);
        }

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
            printf("(%3d, %3d) ", temp->vertexIndex, temp->bondID);
            temp = temp->next;
        }
    }
}

void deleteGraph(AdjacencyList* graph){
    if(graph == NULL)
        return;

    for(int i = 0; i < graph->nodeCount; i++){
        deleteLL(graph->arr[i].next);
    }

    free(graph->arr);
    free(graph);
}

// functions for linked list
// creates a new node of Linked List
SllNode* createSllNode(int vi, int bondID)
{
    SllNode* newNode = (SllNode*)malloc(sizeof(SllNode));
    newNode->vertexIndex = vi;
    newNode->bondID = bondID;
    newNode->next = NULL;
    return newNode;
}

// Function to add new node to the beginning of the Linked List
SllNode* addNodeBegin(SllNode* head, int vi, int bondID)
{
    SllNode* newNode = createSllNode(vi, bondID);
    newNode->next = head;
    return newNode;
}

// Deletes a Linked List
void deleteLL(SllNode* head){
    SllNode* temp, *next;
    temp = head;
    while(temp != NULL){
        next = temp->next;
        free(temp);
        temp = next;
    }
}