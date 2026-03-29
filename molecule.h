#pragma once

#include<stdio.h>
#include<string.h>
#include<stdbool.h>
#include"adjlist.h"

#define COVALENT_BOND 0
#define IONIC_BOND 1

#define SINGLE_BOND 1
#define DOUBLE_BOND 2
#define TRIPLE_BOND 3

typedef struct Atom{
    // int index;
    char element[4];
    int charge;
} Atom;

typedef struct Bond{
    // 1 for single bond; 2 for double bond
    int bondOrder;

    // bondType = 0 (for covalent)
    //          = 1 (for ionic)
    int bondType;
} Bond;

typedef struct Molecule{
    AdjacencyList* graph;
    Atom* atoms;
    Bond* bonds;
} Molecule;

Molecule* createNewMolecule(){
    Molecule* newm = (Molecule *)malloc(sizeof(Molecule));
    newm->graph = createNewAdjacencyList(0);

    newm->atoms = (Atom *)malloc(sizeof(Atom) * newm->graph->capacity);

    newm->bonds = (Bond *)malloc(sizeof(Bond) * newm->graph->bondCapacity);

    return newm;
}

void addAtom(Molecule* mol, char atom[], int charge){
    if(mol->graph->capacity == mol->graph->nodeCount){
        Atom* newAtoms = (Atom *)realloc(mol->atoms, (mol->graph->capacity + SPARE) * sizeof(Atom));
        if(newAtoms == NULL){
            printf("\nMemory Allocation Failed!");
            exit(1);
        }
        mol->atoms = newAtoms;
    }
    int index = addVertexAL(mol->graph);
    strncpy(mol->atoms[index].element, atom, 3);
    mol->atoms[index].charge = charge;
    mol->atoms[index].element[3] = '\0';
}

// bondType = 0 for covalent and 1 for ionic
void addBond(Molecule* mol, int ind1, int ind2, int bondOrder, int bondType){
    if(mol->graph->bondCapacity == mol->graph->bondCount){
        mol->graph->bondCapacity += SPARE;

        Bond* newArr = (Bond *)realloc(mol->bonds, sizeof(Bond) * mol->graph->bondCapacity);

        if(newArr == NULL){
            printf("\nMemory Allocation Failed!");
            exit(1);
        }

        mol->bonds = newArr;
    }

    int bondID = addEdgeAL(mol->graph, ind1, ind2);

    if(bondID == -1){
        printf("\nInvalid bond");
        return;
    }

    mol->bonds[bondID].bondOrder = bondOrder;
    mol->bonds[bondID].bondType = bondType;
}

void deleteMolecule(Molecule* mol){
    if(mol == NULL)
        return;
    
    deleteGraph(mol->graph);
    free(mol->atoms);
    free(mol->bonds);
    free(mol);
}

void printMolecule(Molecule* mol){
    bool *isVisited = (bool *)calloc(mol->graph->bondCount, sizeof(bool));
    SllNode* temp;
    int bondID;

    for(int i = 0; i < mol->graph->nodeCount; i++){
        temp = mol->graph->arr[i].next;

        while(temp != NULL){

            bondID = temp->bondID;

            // if bond not already printed
            if(!isVisited[bondID]){
                
                if(mol->bonds[bondID].bondType == COVALENT_BOND){
                    printf("\n%s(%d) %s(%d)", mol->atoms[i].element, i, mol->atoms[temp->vertexIndex].element, temp->vertexIndex);
                    printf("\tCovalent Bond (Bond order = %d)", mol->bonds[bondID].bondOrder);
                }

                // if ionic bond
                else{
                    printf("\n%s(%d)(charge = %d) ", mol->atoms[i].element, i, mol->atoms[i].charge);
                    printf("%s(%d)(charge = %d)", mol->atoms[temp->vertexIndex].element, temp->vertexIndex, mol->atoms[temp->vertexIndex].charge);
                    printf("\tIonic Bond (%d)", mol->bonds[bondID].bondOrder);
                }
                
                // marking the bond as visited
                isVisited[bondID] = true;
            }

            temp = temp->next;
        }
    }

    free(isVisited);
}