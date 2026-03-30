#pragma once

#include<stdio.h>
#include<stdbool.h>
#include"molecule.h"

int dfsHelperForCycle(AdjacencyList* graph, int start, bool* isVisited, int parent);
int doesMoleculeContainCycle(Molecule* mol);

int dfsHelperForCycle(AdjacencyList* graph, int start, bool* isVisited, int parent){
    SllNode* temp = graph->arr[start].next;
    isVisited[start] = true;

    while(temp != NULL){

        if(isVisited[temp->vertexIndex]){
            
            // if neighbour is already visited and is not parent
            if(temp->vertexIndex != parent)
                return 1; // cycle detected
        }
        else{
            if(dfsHelperForCycle(graph, temp->vertexIndex, isVisited, start))
                return 1;
        }

        temp = temp->next;
    }

    return 0;
}

int doesMoleculeContainCycle(Molecule* mol){
    if(mol == NULL || mol->graph->nodeCount <= 1)
        return 0;
    bool* isVisited = (bool *)calloc(mol->graph->nodeCount, sizeof(bool));

    int res = dfsHelperForCycle(mol->graph, 0, isVisited, -1);
    free(isVisited);

    return res;
}