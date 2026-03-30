#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "molecule.h"

#define MAX_MOLECULES 10

typedef struct Reaction {
    Molecule* reactants[MAX_MOLECULES];
    int reactantCount;

    Molecule* products[MAX_MOLECULES];
    int productCount;

    char name[100];
} Reaction;


// Create Reaction
static Reaction* createReaction(const char* name) {
    Reaction* r = (Reaction*) malloc(sizeof(Reaction));

    r->reactantCount = 0;
    r->productCount = 0;

    strcpy(r->name, name);

    return r;
}


// Add Reactant
static void addReactant(Reaction* r, Molecule* m) {
    if(r->reactantCount < MAX_MOLECULES) {
        r->reactants[r->reactantCount++] = m;
    }
}


// Add Product
static void addProduct(Reaction* r, Molecule* m) {
    if(r->productCount < MAX_MOLECULES) {
        r->products[r->productCount++] = m;
    }
}


// Print Reaction (simple format)
static void printReaction(Reaction* r) {
    printf("\n=== Reaction: %s ===\n", r->name);

    printf("Reactants:\n");
    for(int i = 0; i < r->reactantCount; i++) {
        printMolecule(r->reactants[i]);
        printf("\n");
    }

    printf("-----> Products ----->\n");

    for(int i = 0; i < r->productCount; i++) {
        printMolecule(r->products[i]);
        printf("\n");
    }
}


// Delete Reaction
static void deleteReaction(Reaction* r) {
    free(r);
}