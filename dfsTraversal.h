#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "molecule.h"

/*
 * Iterative DFS traversal using a manual stack (array-based).
 * Prints the atom-by-atom traversal order with bond info.
 * This avoids recursion and explicitly shows the stack state.
 */

#define STACK_MAX 100

typedef struct {
    int items[STACK_MAX];
    int top;
} Stack;

void stackInit(Stack* s) {
    s->top = -1;
}

void stackPush(Stack* s, int val) {
    if (s->top < STACK_MAX - 1)
        s->items[++(s->top)] = val;
}

int stackPop(Stack* s) {
    if (s->top == -1) return -1;
    return s->items[(s->top)--];
}

int stackIsEmpty(Stack* s) {
    return s->top == -1;
}

/*
 * Iterative DFS from a starting atom.
 * Prints each atom visited and its bond connections.
 */
void dfsIterative(Molecule* mol, int start) {
    int n = mol->graph->nodeCount;
    bool* visited = (bool*)calloc(n, sizeof(bool));
    int* parent   = (int*)malloc(n * sizeof(int));

    for (int i = 0; i < n; i++) parent[i] = -1;

    Stack s;
    stackInit(&s);
    stackPush(&s, start);

    printf("\n[Iterative DFS] Starting from atom %s(%d)\n", mol->atoms[start].element, start);
    printf("--------------------------------------------\n");

    int step = 1;

    while (!stackIsEmpty(&s)) {
        int curr = stackPop(&s);

        if (visited[curr]) continue;
        visited[curr] = true;

        // Print current atom
        if (parent[curr] == -1)
            printf("Step %2d | Visit %s(%d)  [start]\n", step++, mol->atoms[curr].element, curr);
        else
            printf("Step %2d | Visit %s(%d)  via bond from %s(%d)\n",
                   step++,
                   mol->atoms[curr].element, curr,
                   mol->atoms[parent[curr]].element, parent[curr]);

        // Push all unvisited neighbors
        SllNode* temp = mol->graph->arr[curr].next;
        while (temp != NULL) {
            int nbr = temp->vertexIndex;
            if (!visited[nbr]) {
                stackPush(&s, nbr);
                if (parent[nbr] == -1)
                    parent[nbr] = curr;  // record how we found it
            }
            temp = temp->next;
        }
    }

    printf("--------------------------------------------\n");
    printf("Total atoms visited: %d / %d\n", step - 1, n);

    free(visited);
    free(parent);
}

/*
 * Prints the molecular formula by counting each element.
 * E.g. Benzene → C6H6, Ethanol → C2H6O
 */
void printMolecularFormula(Molecule* mol) {
    // Supported elements in order (standard chemical formula order)
    const char* order[] = {"C", "H", "N", "O", "S", "P", "F", "Cl", "Br", "Na", "K", "Ca"};
    int nOrder = 12;

    printf("\n[Molecular Formula]\n");
    printf("  ");

    for (int e = 0; e < nOrder; e++) {
        int count = 0;
        for (int i = 0; i < mol->graph->nodeCount; i++) {
            if (strcmp(mol->atoms[i].element, order[e]) == 0)
                count++;
        }
        if (count == 1)
            printf("%s", order[e]);
        else if (count > 1)
            printf("%s%d", order[e], count);
    }
    printf("\n");
}

/*
 * Checks if the molecule graph is connected (no isolated fragments).
 * Uses BFS from atom 0 and counts visited nodes.
 */
void checkConnectivity(Molecule* mol) {
    int n = mol->graph->nodeCount;
    bool* visited = (bool*)calloc(n, sizeof(bool));

    // Simple BFS using an array queue
    int* queue = (int*)malloc(n * sizeof(int));
    int front = 0, rear = 0;

    visited[0] = true;
    queue[rear++] = 0;

    while (front < rear) {
        int curr = queue[front++];
        SllNode* temp = mol->graph->arr[curr].next;
        while (temp != NULL) {
            if (!visited[temp->vertexIndex]) {
                visited[temp->vertexIndex] = true;
                queue[rear++] = temp->vertexIndex;
            }
            temp = temp->next;
        }
    }

    int visitedCount = 0;
    for (int i = 0; i < n; i++)
        if (visited[i]) visitedCount++;

    printf("\n[Connectivity Check]\n");
    if (visitedCount == n)
        printf("  Molecule is CONNECTED (%d atoms, all reachable)\n", n);
    else
        printf("  Molecule is DISCONNECTED (%d / %d atoms reachable from atom 0)\n", visitedCount, n);

    free(visited);
    free(queue);
}

/*
 * Prints degree (number of bonds) of each non-hydrogen atom.
 * Useful for understanding valency.
 */
void printAtomDegrees(Molecule* mol) {
    printf("\n[Atom Degrees (bond count per atom)]\n");
    for (int i = 0; i < mol->graph->nodeCount; i++) {
        if (strcmp(mol->atoms[i].element, "H") == 0) continue;

        int degree = 0;
        SllNode* temp = mol->graph->arr[i].next;
        while (temp != NULL) { degree++; temp = temp->next; }

        printf("  %s(%d) : %d bond%s\n",
               mol->atoms[i].element, i,
               degree, degree == 1 ? "" : "s");
    }
}
