#pragma once
#include <stdio.h>
#include <stdlib.h>
#include "molecule.h"
#include "FGdetection.h"
// target is taking two corresponding array by 1-0 logic print name of fg in priorty which are present
void priority(Molecule *mol){
// priority order: COOH > CHO > OH > CONH2 > CO
    int fg[5] = {0};

    fg[0] = isCarboxylic(mol);
    fg[1] = isAldehyde(mol);
    fg[2] = isAlcohol(mol);
    fg[3] = isAmide(mol);     
    fg[4] = isKetone(mol);

    // corresponding names
    char* names[5] = {"COOH", "CHO", "OH", "CONH2", "CO"};

    // check if no functional group
    int sum = 0;
    for(int i = 0; i < 5; i++){
        sum += fg[i];
    }
    if(sum == 0){
        printf("No functional group detected\n"); // no functional group is present
        return;
    }

    // print in priority order
    printf("Functional groups present in priority order : ");
    for(int i = 0; i < 5; i++){
        if(fg[i] == 1){
            printf("%s ", names[i]);
        }
    }
    printf("\n");
}
