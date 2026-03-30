#pragma once

#ifndef LONGEST_CHAIN_H
#define LONGEST_CHAIN_H

#include "molecule.h"
#include <string.h>

// Helper DFS function
void dfsLongest(Molecule* mol, int curr, bool* visited, int length, int* maxLen) {
    
    visited[curr] = true;

    // update max length
    if(length > *maxLen) {
        *maxLen = length;
    }

    SllNode* temp = mol->graph->arr[curr].next;

    while(temp != NULL) {
        int next = temp->vertexIndex;

        // visit only Carbon atoms and unvisited
        if(!visited[next] && strcmp(mol->atoms[next].element, "C") == 0) {
            dfsLongest(mol, next, visited, length + 1, maxLen);
        }

        temp = temp->next;
    }

    // backtrack
    visited[curr] = false;
}


// Main function to find longest carbon chain
int findLongestCarbonChain(Molecule* mol) {

    if(mol == NULL) return 0;

    int maxLen = 0;

    bool* visited = (bool*)calloc(mol->graph->nodeCount, sizeof(bool));

    for(int i = 0; i < mol->graph->nodeCount; i++) {

        // start only from Carbon atoms
        if(strcmp(mol->atoms[i].element, "C") == 0) {
            dfsLongest(mol, i, visited, 1, &maxLen);
        }
    }

    free(visited);
    return maxLen;
}

#endif