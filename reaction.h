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

static void printReaction(Reaction* r) {
    printf("\n=== Reaction: %s ===\n", r->name);

    // Reactants
    for(int i = 0; i < r->reactantCount; i++) {
        printf("[R%d]\n", i+1);
        printMolecule(r->reactants[i]);
        if(i != r->reactantCount - 1){
            printf(" + ");
        }
    }

    printf("\n--------->\n");

    // Products
    for(int i = 0; i < r->productCount; i++) {
        printf("[P%d]\n", i+1);
        printMolecule(r->products[i]);
        if(i != r->productCount - 1)
            printf(" + ");
    }

    printf("\n========================\n");
}

static void deleteReaction(Reaction* r) {
    if(r == NULL) return;

    // free reactants
    for(int i = 0; i < r->reactantCount; i++) {
        deleteMolecule(r->reactants[i]);
    }

    // free products
    for(int i = 0; i < r->productCount; i++) {
        deleteMolecule(r->products[i]);
    }

    free(r);
}

