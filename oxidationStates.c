#include<stdio.h>
#include<stdbool.h>
#include"molecule.h"

int getENOrderIndex(char atom[]){
    // EN series (more en element first)
    int size = 11;
    char* enSeries[] = {"F", "O", "Cl", "N", "Br", "S", "C", "H", "P", "Na", "K"};

    for(int i = 0; i < size; i++){
        if(!strcmp(atom, enSeries[i])){
            return i;
        }
    }

    // element not found
    return -1;
}

int* calculateOxidationStates(Molecule* mol){
    if(mol == NULL) return NULL;
    bool *isVisited = (bool *)calloc(mol->graph->bondCount, sizeof(bool));

    int* oxidationStates = (int *)malloc(mol->graph->nodeCount * sizeof(int));

    for(int i=0; i < mol->graph->nodeCount; i++){

        // for ionic bonds we just add the charge on it to the oxidation state
        // otherwise the charge is just 0
        oxidationStates[i] = mol->atoms[i].charge;
    }

    SllNode* temp;
    int bondID, en1, en2;

    for(int i = 0; i < mol->graph->nodeCount; i++){


        temp = mol->graph->arr[i].next;

        while(temp != NULL){
            bondID = temp->bondID;

            if( !isVisited[bondID] ){

                if(mol->bonds[bondID].bondType == COVALENT_BOND){

                    // Index in the EN Order for atom i
                    en1 = getENOrderIndex(mol->atoms[i].element);

                    // Index in the EN Order for atom temp->vertexIndex (neighbour of i)
                    en2 = getENOrderIndex(mol->atoms[temp->vertexIndex].element);

                    if(en1 == -1 || en2 == -1){
                        int index = (en1 == -1) ? i : temp->vertexIndex;
                        printf("\n\nThe element %s(%d) is not present in the EN Series. Please add it first", mol->atoms[index].element, index);
                        free(oxidationStates);
                        free(isVisited);
                        return NULL;
                    }

                    // if temp->vertexIndex more EN than i
                    if(en1 > en2){
                        oxidationStates[temp->vertexIndex] -= mol->bonds[bondID].bondOrder;
                        oxidationStates[i] += mol->bonds[bondID].bondOrder;
                    }

                    // if i more EN than temp->vertexIndex
                    else if(en1 < en2){
                        oxidationStates[i] -= mol->bonds[bondID].bondOrder;
                        oxidationStates[temp->vertexIndex] += mol->bonds[bondID].bondOrder;
                    }

                    // don't do anything if same elements
                }

            }

            // no need to change the oxidation states for ionic bonds
            // as we have already changed the oxidation state due to them
            // during initialization of oxidationStates array

            // marking the bond as visited
            isVisited[bondID] = true;

            temp = temp->next;
        }
    }

    free(isVisited);
    return oxidationStates;
}