#include<stdio.h>
#include<string.h>
#include"0801cs241079_adjlist.h"
// #define INITIAL_SIZE 10

typedef struct Atom{
    // int index;
    char element[4];
} Atom;

typedef struct Molecule{
    AdjacencyList* graph;
    Atom* atoms;
} Molecule;


Molecule* createNewMolecule(){
    Molecule* newm = (Molecule *)malloc(sizeof(Molecule));
    newm->graph = createNewAdjacencyList(0, 1);
    newm->atoms = (Atom *)malloc(sizeof(Atom) * newm->graph->capacity);
    return newm;
}

void addAtom(Molecule* mol, char atom[]){
    if(mol->graph->capacity == mol->graph->nodeCount){
        Atom* newAtoms = (Atom *)realloc(mol->atoms, (mol->graph->capacity + SPARE) * sizeof(Atom));
        mol->atoms = newAtoms;
    }
    int index = addVertexAL(mol->graph);
    strcpy(mol->atoms[index].element, atom);
}

void addBond(Molecule* mol, int ind1, int ind2, int bondOrder){
    addEdgeAL(mol->graph, ind1, ind2, bondOrder);
}