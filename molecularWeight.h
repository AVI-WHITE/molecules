#pragma once

#include<stdio.h>
#include<string.h>
#include"molecule.h"

struct AtomWeight{
    char element[4];
    double atomicWt;
};

double getAtomicWt(struct AtomWeight array[], int sizeArray, char* element){
    for(int i = 0; i < sizeArray; i++){
        if(!strcmp(array[i].element, element)){
            return array[i].atomicWt;
        }
    }
    return -1;
}

double molecularWeight(Molecule* mol){

    // array of atomic weights
    // limited elements added for now
    struct AtomWeight weightArray[] = {
        {"H", 1.008}, {"C", 12.011}, {"N", 14.007}, {"O", 15.999},
        {"F", 18.998}, {"Na", 22.99}, {"Mg", 24.305}, {"Al", 26.98}, 
        {"P", 30.974}, {"S", 32.06}, {"Cl", 35.5}, {"K", 39.098}, 
        {"Ca", 40.08}, {"Br", 79.904}, {"I", 126.90},
    };

    int size = sizeof(weightArray) / sizeof(struct AtomWeight);

    double moleculeWt = 0, atomWt;
    for(int i=0; i < mol->graph->nodeCount; i++){
        atomWt = getAtomicWt(weightArray, size, mol->atoms[i].element);
        if(atomWt == -1)
            return -1;
        moleculeWt += atomWt;
    }

    return moleculeWt;
}