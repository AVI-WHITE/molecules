#pragma once
#include <stdio.h>
#include <stdlib.h>
#include "molecule.h"
#include "cycleDetection.h"

// cycle detection for aromatic compounds
// dfsHelperForCycle and doesMoleculeContainCycle defined in cycleDetection.h

// check aromatic compounds here
bool checkAromatic(Molecule* mol){
    if(!doesMoleculeContainCycle(mol)) {
        printf("The compounds does not contain cycle");
        return false; // cycle detection condition checking here
        }
    bool *isVisited = (bool *)calloc(mol->graph->bondCount, sizeof(bool));
    int bondArr[mol->graph->bondCount];
    int idx = 0;

    SllNode* temp;

    // Build bond array (only C-C bonds are taken as they are the main part of carbon chain)
    for(int i = 0; i < mol->graph->nodeCount; i++){

        // skip non-carbon atoms
        if(strcmp(mol->atoms[i].element, "C") != 0) // if C is not that element skip the full iteration
            continue;

        temp = mol->graph->arr[i].next;

        while(temp != NULL){
            int j = temp->vertexIndex;

            // considering only C-C bonds
            if(strcmp(mol->atoms[j].element, "C") != 0){
                temp = temp->next;
                continue;
            }
            int bondID = temp->bondID;

            if(!isVisited[bondID]){ // logic is 1 for double bond , 0 for single bond
                if(mol->bonds[bondID].bondOrder == DOUBLE_BOND)
                    bondArr[idx++] = 1;
                else
                    bondArr[idx++] = 0;

                isVisited[bondID] = true;
            }

            temp = temp->next;
        }
    }

    free(isVisited);

    int n = idx;

    //Checking alternating pattern for conjugation of double bonds
    for(int i = 0; i < n-1; i++){
        if(bondArr[i] == bondArr[i+1]){
            printf("there must be no consecutive single or double bond\n");
            return false;
        }
    }

    if(bondArr[0] == bondArr[n-1]){ // if connection of first and last element or we can say where cycle is connected
        printf("there must be no consecutive single or double bond between first and last element\n");
        return false;
    }

    // Checking planarity for Carbon atoms
    for(int i = 0; i < mol->graph->nodeCount; i++){

        if(strcmp(mol->atoms[i].element, "C") != 0)
            continue;

        int hasDouble = 0;
        temp = mol->graph->arr[i].next;

        while(temp != NULL){
            int bondID = temp->bondID;

            if(mol->bonds[bondID].bondOrder == DOUBLE_BOND){ // if we get single double bond in adjacency list
                hasDouble = 1; // of any carbon then it will not be sp3
                break;
            }

            temp = temp->next;
        }

        if(!hasDouble){
            printf("carbon atom has only single bonds means sp3 and will be non planar\n");
            return false;
        }
    }

    // Validating Hückel Rule (π electrons)
    int piee = 0;
    for(int i = 0; i < n; i++){
        if(bondArr[i] == 1) piee += 2; // counting double bonds 
    }

    if((piee - 2) % 4 != 0){ // checking by huckel rule
        printf("the molecule does not follows huckel rule it has not 4n+2 pi electrons\n");
        return false;
    }

    return true;
}

