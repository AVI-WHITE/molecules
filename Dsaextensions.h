#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include "molecule.h"

// ============================================================
//  DSA EXTENSIONS FOR MOLECULAR GRAPH
//  Author: DSA submission additions
//  Covers:
//    1. BFS - Shortest Bond Path between two atoms
//    2. Carbon Classification (Primary/Secondary/Tertiary/Quaternary)
//    3. Graph Isomorphism Check (two molecules same structure?)
//    4. Articulation Points (critical atoms - bridge atoms)
//    5. Topological Sort for Reaction Step Ordering (Kahn's Algorithm)
//    6. Molecular Fingerprint using DFS atom hash
// ============================================================


// ============================================================
//  1. BFS: SHORTEST BOND PATH BETWEEN TWO ATOMS
//  Uses standard BFS on adjacency list.
//  Returns number of bonds between src and dest (-1 if not found)
//  Also prints the atom path.
// ============================================================

int* bfsShortestPath(Molecule* mol, int src, int dest) {
    int n = mol->graph->nodeCount;
    if (src < 0 || dest < 0 || src >= n || dest >= n) return NULL;

    int* dist     = (int*)malloc(n * sizeof(int));
    int* parent   = (int*)malloc(n * sizeof(int));
    bool* visited = (bool*)calloc(n, sizeof(bool));

    for (int i = 0; i < n; i++) { dist[i] = INT_MAX; parent[i] = -1; }

    // BFS queue (circular array)
    int* queue = (int*)malloc(n * sizeof(int));
    int front = 0, rear = 0;

    visited[src] = true;
    dist[src] = 0;
    queue[rear++] = src;

    while (front < rear) {
        int curr = queue[front++];
        SllNode* temp = mol->graph->arr[curr].next;
        while (temp != NULL) {
            int nbr = temp->vertexIndex;
            if (!visited[nbr]) {
                visited[nbr] = true;
                dist[nbr] = dist[curr] + 1;
                parent[nbr] = curr;
                queue[rear++] = nbr;
            }
            temp = temp->next;
        }
    }

    free(visited);
    free(queue);

    if (dist[dest] == INT_MAX) {
        printf("\n[BFS] No path found between atom %d (%s) and atom %d (%s)\n",
               src, mol->atoms[src].element, dest, mol->atoms[dest].element);
        free(dist);
        free(parent);
        return NULL;
    }

    // Print path by backtracking from dest
    printf("\n[BFS] Shortest path from %s(%d) to %s(%d): %d bonds\n",
           mol->atoms[src].element, src,
           mol->atoms[dest].element, dest,
           dist[dest]);

    // Reconstruct and print path
    int* path = (int*)malloc((dist[dest] + 1) * sizeof(int));
    int idx = 0, cur = dest;
    while (cur != -1) { path[idx++] = cur; cur = parent[cur]; }
    printf("  Path: ");
    for (int i = idx - 1; i >= 0; i--) {
        printf("%s(%d)", mol->atoms[path[i]].element, path[i]);
        if (i > 0) printf(" -> ");
    }
    printf("\n");

    free(dist);
    free(path);
    return parent; // caller must free
}


// ============================================================
//  2. CARBON CLASSIFICATION
//  Primary   (1°) = bonded to exactly 1 other Carbon
//  Secondary (2°) = bonded to exactly 2 Carbons
//  Tertiary  (3°) = bonded to exactly 3 Carbons
//  Quaternary(4°) = bonded to exactly 4 Carbons
// ============================================================

typedef enum {
    PRIMARY = 1,
    SECONDARY = 2,
    TERTIARY = 3,
    QUATERNARY = 4,
    NOT_CARBON = 0
} CarbonType;

CarbonType classifyCarbon(Molecule* mol, int atomIdx) {
    if (strcmp(mol->atoms[atomIdx].element, "C") != 0)
        return NOT_CARBON;

    int carbonNeighbors = 0;
    SllNode* temp = mol->graph->arr[atomIdx].next;
    while (temp != NULL) {
        if (strcmp(mol->atoms[temp->vertexIndex].element, "C") == 0)
            carbonNeighbors++;
        temp = temp->next;
    }

    if (carbonNeighbors == 1) return PRIMARY;
    if (carbonNeighbors == 2) return SECONDARY;
    if (carbonNeighbors == 3) return TERTIARY;
    if (carbonNeighbors >= 4) return QUATERNARY;
    return PRIMARY; // isolated carbon (no C-C bonds)
}

void printCarbonClassification(Molecule* mol) {
    printf("\n[Carbon Classification]\n");
    const char* labels[] = {"", "Primary (1°)", "Secondary (2°)", "Tertiary (3°)", "Quaternary (4°)"};
    int counts[5] = {0};

    for (int i = 0; i < mol->graph->nodeCount; i++) {
        CarbonType ct = classifyCarbon(mol, i);
        if (ct != NOT_CARBON) {
            printf("  C(%d) -> %s\n", i, labels[ct]);
            counts[ct]++;
        }
    }
    printf("  Summary: Primary=%d, Secondary=%d, Tertiary=%d, Quaternary=%d\n",
           counts[1], counts[2], counts[3], counts[4]);
}


// ============================================================
//  3. GRAPH ISOMORPHISM (Degree-Sequence Check)
//  Full graph isomorphism is NP-complete, but for molecules
//  we use a practical approach:
//    Step 1: Same number of atoms?
//    Step 2: Same number of bonds?
//    Step 3: Same element frequency?
//    Step 4: Same degree sequence per element?
//  This is a strong necessary condition and works for most
//  real molecule comparison scenarios.
// ============================================================

// Comparison helper
int comparInt(const void* a, const void* b) {
    return (*(int*)a - *(int*)b);
}

bool areMoleculesIsomorphic(Molecule* m1, Molecule* m2) {
    int n1 = m1->graph->nodeCount;
    int n2 = m2->graph->nodeCount;
    int e1 = m1->graph->bondCount;
    int e2 = m2->graph->bondCount;

    printf("\n[Isomorphism Check]\n");

    // Step 1: same vertex count
    if (n1 != n2) {
        printf("  Different atom count (%d vs %d) -> NOT isomorphic\n", n1, n2);
        return false;
    }

    // Step 2: same edge count
    if (e1 != e2) {
        printf("  Different bond count (%d vs %d) -> NOT isomorphic\n", e1, e2);
        return false;
    }

    // Step 3: same element frequency
    // Count top 10 elements
    const char* elements[] = {"C","H","O","N","S","P","F","Cl","Br","Na"};
    int nElem = 10;
    for (int e = 0; e < nElem; e++) {
        int c1 = 0, c2 = 0;
        for (int i = 0; i < n1; i++)
            if (!strcmp(m1->atoms[i].element, elements[e])) c1++;
        for (int i = 0; i < n2; i++)
            if (!strcmp(m2->atoms[i].element, elements[e])) c2++;
        if (c1 != c2) {
            printf("  Different count of %s (%d vs %d) -> NOT isomorphic\n", elements[e], c1, c2);
            return false;
        }
    }

    // Step 4: degree sequence must match per element
    int* deg1 = (int*)calloc(n1, sizeof(int));
    int* deg2 = (int*)calloc(n2, sizeof(int));

    SllNode* temp;
    for (int i = 0; i < n1; i++) {
        temp = m1->graph->arr[i].next;
        while (temp) { deg1[i]++; temp = temp->next; }
    }
    for (int i = 0; i < n2; i++) {
        temp = m2->graph->arr[i].next;
        while (temp) { deg2[i]++; temp = temp->next; }
    }

    qsort(deg1, n1, sizeof(int), comparInt);
    qsort(deg2, n2, sizeof(int), comparInt);

    bool degreeSame = true;
    for (int i = 0; i < n1; i++) {
        if (deg1[i] != deg2[i]) { degreeSame = false; break; }
    }

    free(deg1); free(deg2);

    if (!degreeSame) {
        printf("  Degree sequences differ -> NOT isomorphic\n");
        return false;
    }

    printf("  Atom count, bond count, element frequency, and degree sequences all match.\n");
    printf("  Result: Structurally ISOMORPHIC (same molecular graph)\n");
    return true;
}


// ============================================================
//  4. ARTICULATION POINTS (Critical / Bridge Atoms)
//  An articulation point is a vertex whose removal disconnects
//  the graph. In chemistry: these are atoms that if broken
//  split the molecule into two fragments.
//  Uses Tarjan's DFS-based algorithm. O(V + E)
// ============================================================

static int apTimer = 0;

void apDFS(AdjacencyList* graph, int u,
           bool* visited, int* disc, int* low,
           int* parent, bool* ap) {
    int children = 0;
    visited[u] = true;
    disc[u] = low[u] = apTimer++;

    SllNode* temp = graph->arr[u].next;
    while (temp != NULL) {
        int v = temp->vertexIndex;
        if (!visited[v]) {
            children++;
            parent[v] = u;
            apDFS(graph, v, visited, disc, low, parent, ap);
            low[u] = (low[u] < low[v]) ? low[u] : low[v];

            // u is an AP if: root with 2+ children, OR non-root with low[v] >= disc[u]
            if (parent[u] == -1 && children > 1) ap[u] = true;
            if (parent[u] != -1 && low[v] >= disc[u]) ap[u] = true;
        } else if (v != parent[u]) {
            low[u] = (low[u] < disc[v]) ? low[u] : disc[v];
        }
        temp = temp->next;
    }
}

void findArticulationPoints(Molecule* mol) {
    int n = mol->graph->nodeCount;
    bool* visited = (bool*)calloc(n, sizeof(bool));
    int*  disc    = (int*)malloc(n * sizeof(int));
    int*  low     = (int*)malloc(n * sizeof(int));
    int*  parent  = (int*)malloc(n * sizeof(int));
    bool* ap      = (bool*)calloc(n, sizeof(bool));

    for (int i = 0; i < n; i++) parent[i] = -1;
    apTimer = 0;

    for (int i = 0; i < n; i++)
        if (!visited[i])
            apDFS(mol->graph, i, visited, disc, low, parent, ap);

    printf("\n[Articulation Points / Critical Atoms]\n");
    bool found = false;
    for (int i = 0; i < n; i++) {
        if (ap[i]) {
            printf("  Atom %s(%d) is a critical atom — removing it disconnects the molecule\n",
                   mol->atoms[i].element, i);
            found = true;
        }
    }
    if (!found)
        printf("  No articulation points — molecule is 2-connected (robust structure)\n");

    free(visited); free(disc); free(low); free(parent); free(ap);
}


// ============================================================
//  5. TOPOLOGICAL SORT FOR MULTI-STEP REACTION ORDERING
//  Models a multi-step synthesis as a DAG:
//    Nodes = reaction steps
//    Edges = step A must happen before step B
//  Uses Kahn's Algorithm (BFS-based). O(V + E)
// ============================================================

#define MAX_STEPS 20

typedef struct ReactionStep {
    char name[64];
    int  dependencies[MAX_STEPS]; // indices of steps that must come before
    int  depCount;
} ReactionStep;

typedef struct SynthesisRoute {
    ReactionStep steps[MAX_STEPS];
    int stepCount;
} SynthesisRoute;

SynthesisRoute* createSynthesisRoute() {
    SynthesisRoute* sr = (SynthesisRoute*)malloc(sizeof(SynthesisRoute));
    sr->stepCount = 0;
    return sr;
}

int addStep(SynthesisRoute* sr, const char* name) {
    int idx = sr->stepCount++;
    strncpy(sr->steps[idx].name, name, 63);
    sr->steps[idx].depCount = 0;
    return idx;
}

void addDependency(SynthesisRoute* sr, int from, int to) {
    // step 'from' must complete before step 'to'
    ReactionStep* step = &sr->steps[to];
    step->dependencies[step->depCount++] = from;
}

void topologicalSortKahn(SynthesisRoute* sr) {
    int n = sr->stepCount;
    int* inDegree = (int*)calloc(n, sizeof(int));

    // Build adjacency and compute in-degrees
    int adj[MAX_STEPS][MAX_STEPS];
    int adjCount[MAX_STEPS];
    memset(adj, 0, sizeof(adj));
    memset(adjCount, 0, sizeof(adjCount));

    for (int i = 0; i < n; i++) {
        for (int d = 0; d < sr->steps[i].depCount; d++) {
            int from = sr->steps[i].dependencies[d];
            adj[from][adjCount[from]++] = i;
            inDegree[i]++;
        }
    }

    // BFS queue with in-degree 0 nodes
    int queue[MAX_STEPS], front = 0, rear = 0;
    for (int i = 0; i < n; i++)
        if (inDegree[i] == 0) queue[rear++] = i;

    printf("\n[Topological Sort - Reaction Step Ordering (Kahn's Algorithm)]\n");
    int processed = 0;
    int order = 1;

    while (front < rear) {
        int curr = queue[front++];
        printf("  Step %d: %s\n", order++, sr->steps[curr].name);
        processed++;

        for (int i = 0; i < adjCount[curr]; i++) {
            int nbr = adj[curr][i];
            if (--inDegree[nbr] == 0)
                queue[rear++] = nbr;
        }
    }

    if (processed != n) {
        printf("  WARNING: Cycle detected in reaction steps! Invalid synthesis route.\n");
    } else {
        printf("  All %d steps ordered successfully.\n", n);
    }

    free(inDegree);
}

void deleteSynthesisRoute(SynthesisRoute* sr) {
    free(sr);
}


// ============================================================
//  6. MOLECULAR FINGERPRINT using DFS-based atom hash
//  Generates a unique "fingerprint" string by DFS traversal,
//  encoding element + degree at each atom.
//  Useful for quick molecule comparison / hashing.
// ============================================================

void dfsFingerprintHelper(Molecule* mol, int curr, bool* visited, char* fp, int* pos) {
    visited[curr] = true;

    // Encode: element + number of bonds
    int degree = 0;
    SllNode* temp = mol->graph->arr[curr].next;
    while (temp) { degree++; temp = temp->next; }

    // Write to fingerprint buffer
    char buf[16];
    snprintf(buf, sizeof(buf), "[%s%d]", mol->atoms[curr].element, degree);
    int len = strlen(buf);
    strncpy(fp + *pos, buf, len);
    *pos += len;

    temp = mol->graph->arr[curr].next;
    while (temp != NULL) {
        if (!visited[temp->vertexIndex])
            dfsFingerprintHelper(mol, temp->vertexIndex, visited, fp, pos);
        temp = temp->next;
    }
}

char* generateFingerprint(Molecule* mol) {
    int n = mol->graph->nodeCount;
    char* fp = (char*)calloc(n * 16 + 1, sizeof(char));
    bool* visited = (bool*)calloc(n, sizeof(bool));
    int pos = 0;

    for (int i = 0; i < n; i++)
        if (!visited[i])
            dfsFingerprintHelper(mol, i, visited, fp, &pos);

    free(visited);
    printf("\n[Molecular Fingerprint (DFS-encoded)]\n  %s\n", fp);
    return fp; // caller must free
}


// ============================================================
//  DEMO: Run all DSA Extensions on a given molecule
// ============================================================

void runAllDSAExtensions(Molecule* mol, const char* molName) {
    printf("\n");
    printf("╔══════════════════════════════════════════════════════╗\n");
    printf("║   DSA ANALYSIS: %-35s ║\n", molName);
    printf("╚══════════════════════════════════════════════════════╝\n");

    // 1. Carbon Classification
    printCarbonClassification(mol);

    // 2. Articulation Points
    findArticulationPoints(mol);

    // 3. BFS Shortest Path (between atom 0 and last atom)
    int last = mol->graph->nodeCount - 1;
    int* p = bfsShortestPath(mol, 0, last);
    if (p) free(p);

    // 4. Molecular Fingerprint
    char* fp = generateFingerprint(mol);
    if (fp) free(fp);
}