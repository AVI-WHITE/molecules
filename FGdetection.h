#pragma once
#include <stdio.h>
#include <stdlib.h>
#include "molecule.h"
// detection of 5 main functional group

// alcohol detection
bool isAlcohol(Molecule* mol){
    SllNode* temp;

    for(int i = 0; i < mol->graph->nodeCount; i++){

        if(strcmp(mol->atoms[i].element, "O") != 0) // for an atom O CONTAINed in adjacency list checking
            continue;
        
        int hasH = 0, hasC = 0;
        int carbonIndex = -1;

        temp = mol->graph->arr[i].next;

        while(temp != NULL){
            int j = temp->vertexIndex;
            // if O atom adjaceny list has H and C has its connected atom then it is alcohol

            if(strcmp(mol->atoms[j].element, "H") == 0 && mol->bonds[temp->bondID].bondOrder == SINGLE_BOND){
                hasH = 1;
            }
            if(strcmp(mol->atoms[j].element, "C") == 0 && mol->bonds[temp->bondID].bondOrder == SINGLE_BOND){
                hasC = 1;
                carbonIndex = j; // storing carbonindex to later differentiate bw carboxylix OH and alcohol
            }

            temp = temp->next;
        }

        if(hasH && hasC){

            //check that this carbon is NOT part of C=O (COOH case)
            int hasCarbonyl = 0;
            temp = mol->graph->arr[carbonIndex].next;

            while(temp != NULL){
                int k = temp->vertexIndex;
                  // i f same carbon which contained OH group has c=o in it then it will become carboxylic acid
                if(strcmp(mol->atoms[k].element, "O") == 0 && mol->bonds[temp->bondID].bondOrder == DOUBLE_BOND){
                    hasCarbonyl = 1; // for double bond O
                    break;
                }

                temp = temp->next;
            }

            // alcohol only if no C=O on that carbon
            if(!hasCarbonyl)
                return true;
        }
    }

    return false;
}


bool isAldehyde(Molecule* mol){
    SllNode* temp;
    for(int i = 0; i < mol->graph->nodeCount; i++){

        if(strcmp(mol->atoms[i].element, "C") != 0) // process only Carbon
            continue;

        int hasOdouble = 0, hasH = 0;
        temp = mol->graph->arr[i].next;

        while(temp != NULL){
            int j = temp->vertexIndex;

            // check C=O that for the carbon only if it has one c=o and h connected to it then it must be aldehyde
            if(strcmp(mol->atoms[j].element, "O") == 0 && mol->bonds[temp->bondID].bondOrder == DOUBLE_BOND){
                hasOdouble = 1;
            }
            // check C-H
            if(strcmp(mol->atoms[j].element, "H") == 0 && mol->bonds[temp->bondID].bondOrder == SINGLE_BOND){
                hasH = 1;
            }

            temp = temp->next;
        }

        if(hasOdouble == 1 && hasH == 1)
            return true;
    }

    return false;
}

bool isCarboxylic(Molecule* mol){
    SllNode* temp;

    for(int i = 0; i < mol->graph->nodeCount; i++){

        if(strcmp(mol->atoms[i].element, "C") != 0)
            continue;

        int hasOdouble = 0, hasOH = 0;
        temp = mol->graph->arr[i].next;

        while(temp != NULL){
            int j = temp->vertexIndex;
            int bondID = temp->bondID;

            // check C=O
            if(strcmp(mol->atoms[j].element, "O") == 0 && mol->bonds[bondID].bondOrder == DOUBLE_BOND){
                hasOdouble = 1;
            }
            // check C-O-H
            if(strcmp(mol->atoms[j].element, "O") == 0 && mol->bonds[bondID].bondOrder == SINGLE_BOND){
                // if carbon has oxygen connected to it then we have to check other connection of oxygen
                SllNode* t2 = mol->graph->arr[j].next;
                while(t2 != NULL){
                    int k = t2->vertexIndex;

                    if(strcmp(mol->atoms[k].element, "H") == 0 && mol->bonds[t2->bondID].bondOrder == SINGLE_BOND){ // if h present
                        hasOH = 1;
                        break;
                    }
                  t2 = t2->next;
                }
            }
            temp = temp->next;
        }
        if(hasOdouble && hasOH) return true; // both has carbonyl and OH then carboxylic acid
    }
    return false;
}


bool isAmide(Molecule* mol){ // we have to check if c is connnected with O and N atom then it will be amide of any degree
    SllNode* temp;

    for(int i = 0; i < mol->graph->nodeCount; i++){

        if(strcmp(mol->atoms[i].element, "C") != 0)
            continue;

        int hasOdouble = 0, hasN = 0;
        temp = mol->graph->arr[i].next;

        while(temp != NULL){
            int j = temp->vertexIndex;
            int bondID = temp->bondID;

            // check C=O
            if(strcmp(mol->atoms[j].element, "O") == 0 && mol->bonds[bondID].bondOrder == DOUBLE_BOND){
                hasOdouble = 1;
            }
            // check C-N-C or C-N-H
            if(strcmp(mol->atoms[j].element, "N") == 0 &&  mol->bonds[bondID].bondOrder == SINGLE_BOND){
                hasN = 1;
            }
            temp = temp->next;
        }
        if(hasOdouble && hasN)
            return true;
    }
    return false;
}


bool isKetone(Molecule* mol){
    SllNode* temp;
    for(int i = 0; i < mol->graph->nodeCount; i++){

        if(strcmp(mol->atoms[i].element, "C") != 0) // process only Carbon
            continue;

        int hasOdouble = 0, hasC = 0,hasDO =0 ;
        temp = mol->graph->arr[i].next;

        while(temp != NULL){
            int j = temp->vertexIndex;

            // check C=O
            if(strcmp(mol->atoms[j].element, "O") == 0 && mol->bonds[temp->bondID].bondOrder == DOUBLE_BOND){ // carbon must have c=o
                hasOdouble = 1;
            }
            // check C-H
            if(strcmp(mol->atoms[j].element, "C") == 0 && mol->bonds[temp->bondID].bondOrder == SINGLE_BOND){ // and other two bonds must be connected with carbon to ensure ketone
                hasC = 2;
            }

            temp = temp->next;
        }

        if(hasOdouble == 1 && hasC == 2)
            return true;
    }

    return false;
}

