/*
 * ============================================================
 *  mainTest.c  —  MOLECULES PROJECT : COMPLETE TEST FILE
 *
 *  Runs ALL features in one place:
 *    PART 1 : Core molecule operations (existing project)
 *             - Benzene, Na2SO4, Isobutane, Naphthalene
 *             - Molecular weight, oxidation states,
 *               longest carbon chain, cycle detection
 *
 *    PART 2 : Functional group detection & priority
 *             - Alcohol, Aldehyde, Ketone, Carboxylic acid, Amide
 *             - Priority ordering (COOH > CHO > OH > CONH2 > CO)
 *
 *    PART 3 : Aromatic detection (Huckel rule)
 *             - Benzene (aromatic), Isobutane (not aromatic)
 *
 *    PART 4 : Reaction representation
 *             - Reactants + Products
 *
 *    PART 5 : NEW DSA ADDITIONS (dfsTraversal.h)
 *             - Iterative DFS with manual stack
 *             - Molecular formula printer
 *             - Connectivity check (BFS)
 *             - Atom degree analysis
 *
 *  Compile : gcc mainTest.c -o mainTest
 *  Run     : ./mainTest
 * ============================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include "molecule.h"
#include "operationsOnMolecules.h"
#include "FGdetection.h"
#include "FGpriority.h"
#include "aromatic.h"
#include "reaction.h"
#include "dfsTraversal.h"

/* ── helper: print a section banner ─────────────────────── */
void banner(const char* title) {
    printf("\n");
    printf("============================================================\n");
    printf("  %s\n", title);
    printf("============================================================\n");
}

void subbanner(const char* title) {
    printf("\n--- %s ---\n", title);
}

/* ============================================================
   PART 1 : CORE OPERATIONS
   ============================================================ */
void part1_core() {
    banner("PART 1 : Core Molecule Operations");

    /* ── Benzene ── */
    subbanner("Benzene (C6H6)");
    Molecule* benzene = createNewMolecule();
    for (int i = 0; i < 6; i++) addAtom(benzene, "C", 0);
    for (int i = 0; i < 6; i++) addAtom(benzene, "H", 0);
    addBond(benzene, 0, 1, DOUBLE_BOND, COVALENT_BOND);
    addBond(benzene, 1, 2, SINGLE_BOND,  COVALENT_BOND);
    addBond(benzene, 2, 3, DOUBLE_BOND, COVALENT_BOND);
    addBond(benzene, 3, 4, SINGLE_BOND,  COVALENT_BOND);
    addBond(benzene, 4, 5, DOUBLE_BOND, COVALENT_BOND);
    addBond(benzene, 5, 0, SINGLE_BOND,  COVALENT_BOND);
    for (int i = 0; i < 6; i++)
        addBond(benzene, i, i + 6, SINGLE_BOND, COVALENT_BOND);

    printMolecule(benzene);
    printf("\nMolecular weight    : %.3f g/mol", molecularWeight(benzene));
    printf("\nLongest carbon chain: %d", findLongestCarbonChain(benzene));
    printf("\nIs cyclic           : %d\n", doesMoleculeContainCycle(benzene));

    /* ── Na2SO4 ── */
    subbanner("Na2SO4");
    Molecule* na2so4 = createNewMolecule();
    addAtom(na2so4, "S",  0);
    addAtom(na2so4, "O",  0);
    addAtom(na2so4, "O",  0);
    addAtom(na2so4, "O", -1);
    addAtom(na2so4, "O", -1);
    addAtom(na2so4, "Na", 1);
    addAtom(na2so4, "Na", 1);
    addBond(na2so4, 0, 1, DOUBLE_BOND, COVALENT_BOND);
    addBond(na2so4, 0, 2, DOUBLE_BOND, COVALENT_BOND);
    addBond(na2so4, 0, 3, SINGLE_BOND, COVALENT_BOND);
    addBond(na2so4, 0, 4, SINGLE_BOND, COVALENT_BOND);
    addBond(na2so4, 3, 5, SINGLE_BOND, IONIC_BOND);
    addBond(na2so4, 4, 6, SINGLE_BOND, IONIC_BOND);

    printMolecule(na2so4);
    printf("\nMolecular weight : %.3f g/mol", molecularWeight(na2so4));
    int* oxNa = calculateOxidationStates(na2so4);
    printOxidationStatesArray(na2so4, oxNa);
    free(oxNa);

    /* ── Isobutane ── */
    subbanner("Isobutane (C4H10)");
    Molecule* isobutane = createNewMolecule();
    addAtom(isobutane, "C", 0); // 0 central
    addAtom(isobutane, "C", 0); // 1
    addAtom(isobutane, "C", 0); // 2
    addAtom(isobutane, "C", 0); // 3
    for (int i = 0; i < 10; i++) addAtom(isobutane, "H", 0);
    addBond(isobutane, 0, 1,  SINGLE_BOND, COVALENT_BOND);
    addBond(isobutane, 0, 2,  SINGLE_BOND, COVALENT_BOND);
    addBond(isobutane, 0, 3,  SINGLE_BOND, COVALENT_BOND);
    addBond(isobutane, 0, 4,  SINGLE_BOND, COVALENT_BOND);
    addBond(isobutane, 1, 5,  SINGLE_BOND, COVALENT_BOND);
    addBond(isobutane, 1, 6,  SINGLE_BOND, COVALENT_BOND);
    addBond(isobutane, 1, 7,  SINGLE_BOND, COVALENT_BOND);
    addBond(isobutane, 2, 8,  SINGLE_BOND, COVALENT_BOND);
    addBond(isobutane, 2, 9,  SINGLE_BOND, COVALENT_BOND);
    addBond(isobutane, 2, 10, SINGLE_BOND, COVALENT_BOND);
    addBond(isobutane, 3, 11, SINGLE_BOND, COVALENT_BOND);
    addBond(isobutane, 3, 12, SINGLE_BOND, COVALENT_BOND);
    addBond(isobutane, 3, 13, SINGLE_BOND, COVALENT_BOND);

    printf("Is cyclic            : %d\n", doesMoleculeContainCycle(isobutane));
    printf("Longest carbon chain : %d\n", findLongestCarbonChain(isobutane));
    printf("Molecular weight     : %.3f g/mol\n", molecularWeight(isobutane));
    int* oxIso = calculateOxidationStates(isobutane);
    printOxidationStatesArray(isobutane, oxIso);
    free(oxIso);

    /* ── Naphthalene ── */
    subbanner("Naphthalene (C10H8)");
    Molecule* naph = createNewMolecule();
    for (int i = 0; i < 10; i++) addAtom(naph, "C", 0);
    for (int i = 0; i < 8;  i++) addAtom(naph, "H", 0);
    addBond(naph, 0, 1, DOUBLE_BOND, COVALENT_BOND);
    addBond(naph, 1, 2, SINGLE_BOND, COVALENT_BOND);
    addBond(naph, 2, 3, DOUBLE_BOND, COVALENT_BOND);
    addBond(naph, 3, 4, SINGLE_BOND, COVALENT_BOND);
    addBond(naph, 4, 9, DOUBLE_BOND, COVALENT_BOND);
    addBond(naph, 9, 0, SINGLE_BOND, COVALENT_BOND);
    addBond(naph, 4, 5, SINGLE_BOND, COVALENT_BOND);
    addBond(naph, 5, 6, DOUBLE_BOND, COVALENT_BOND);
    addBond(naph, 6, 7, SINGLE_BOND, COVALENT_BOND);
    addBond(naph, 7, 8, DOUBLE_BOND, COVALENT_BOND);
    addBond(naph, 8, 9, SINGLE_BOND, COVALENT_BOND);
    addBond(naph, 0, 10, SINGLE_BOND, COVALENT_BOND);
    addBond(naph, 1, 11, SINGLE_BOND, COVALENT_BOND);
    addBond(naph, 2, 12, SINGLE_BOND, COVALENT_BOND);
    addBond(naph, 3, 13, SINGLE_BOND, COVALENT_BOND);
    addBond(naph, 5, 14, SINGLE_BOND, COVALENT_BOND);
    addBond(naph, 6, 15, SINGLE_BOND, COVALENT_BOND);
    addBond(naph, 7, 16, SINGLE_BOND, COVALENT_BOND);
    addBond(naph, 8, 17, SINGLE_BOND, COVALENT_BOND);

    printf("Is cyclic            : %d\n", doesMoleculeContainCycle(naph));
    printf("Longest carbon chain : %d\n", findLongestCarbonChain(naph));
    printf("Molecular weight     : %.3f g/mol\n", molecularWeight(naph));
    int* oxNaph = calculateOxidationStates(naph);
    printOxidationStatesArray(naph, oxNaph);
    free(oxNaph);

    deleteMolecule(benzene);
    deleteMolecule(na2so4);
    deleteMolecule(isobutane);
    deleteMolecule(naph);
}

/* ============================================================
   PART 2 : FUNCTIONAL GROUP DETECTION & PRIORITY
   ============================================================ */
void part2_functional_groups() {
    banner("PART 2 : Functional Group Detection & Priority");

    /* ── Alcohol: Ethanol ── */
    subbanner("Ethanol (Alcohol)");
    Molecule* ethanol = createNewMolecule();
    addAtom(ethanol, "C", 0); addAtom(ethanol, "C", 0);
    addAtom(ethanol, "O", 0);
    addAtom(ethanol, "H", 0); addAtom(ethanol, "H", 0);
    addAtom(ethanol, "H", 0); addAtom(ethanol, "H", 0);
    addAtom(ethanol, "H", 0); addAtom(ethanol, "H", 0);
    addBond(ethanol, 0, 1, SINGLE_BOND, COVALENT_BOND);
    addBond(ethanol, 1, 2, SINGLE_BOND, COVALENT_BOND);
    addBond(ethanol, 2, 8, SINGLE_BOND, COVALENT_BOND);
    addBond(ethanol, 0, 3, SINGLE_BOND, COVALENT_BOND);
    addBond(ethanol, 0, 4, SINGLE_BOND, COVALENT_BOND);
    addBond(ethanol, 0, 5, SINGLE_BOND, COVALENT_BOND);
    addBond(ethanol, 1, 6, SINGLE_BOND, COVALENT_BOND);
    addBond(ethanol, 1, 7, SINGLE_BOND, COVALENT_BOND);
    printf("isAlcohol    : %d\n", isAlcohol(ethanol));
    printf("isAldehyde   : %d\n", isAldehyde(ethanol));
    printf("isKetone     : %d\n", isKetone(ethanol));
    printf("isCarboxylic : %d\n", isCarboxylic(ethanol));
    printf("isAmide      : %d\n", isAmide(ethanol));
    printf("Priority     : "); priority(ethanol);
    deleteMolecule(ethanol);

    /* ── Carboxylic acid: Acetic acid ── */
    subbanner("Acetic Acid (Carboxylic acid)");
    Molecule* acid = createNewMolecule();
    addAtom(acid, "C", 0); // 0 CH3
    addAtom(acid, "C", 0); // 1 COOH carbon
    addAtom(acid, "O", 0); // 2 =O
    addAtom(acid, "O", 0); // 3 OH
    addAtom(acid, "H", 0); addAtom(acid, "H", 0);
    addAtom(acid, "H", 0); addAtom(acid, "H", 0);
    addBond(acid, 0, 1, SINGLE_BOND, COVALENT_BOND);
    addBond(acid, 1, 2, DOUBLE_BOND, COVALENT_BOND);
    addBond(acid, 1, 3, SINGLE_BOND, COVALENT_BOND);
    addBond(acid, 3, 7, SINGLE_BOND, COVALENT_BOND);
    addBond(acid, 0, 4, SINGLE_BOND, COVALENT_BOND);
    addBond(acid, 0, 5, SINGLE_BOND, COVALENT_BOND);
    addBond(acid, 0, 6, SINGLE_BOND, COVALENT_BOND);
    printf("isCarboxylic : %d\n", isCarboxylic(acid));
    printf("Priority     : "); priority(acid);
    deleteMolecule(acid);

    /* ── Aldehyde: Butanal ── */
    subbanner("Butanal (Aldehyde)");
    Molecule* aldehyde = createNewMolecule();
    addAtom(aldehyde, "C", 0); // 0
    addAtom(aldehyde, "C", 0); // 1
    addAtom(aldehyde, "C", 0); // 2
    addAtom(aldehyde, "C", 0); // 3 CHO carbon
    addAtom(aldehyde, "O", 0); // 4 =O
    addAtom(aldehyde, "H", 0); addAtom(aldehyde, "H", 0);
    addAtom(aldehyde, "H", 0); addAtom(aldehyde, "H", 0);
    addAtom(aldehyde, "H", 0); addAtom(aldehyde, "H", 0);
    addAtom(aldehyde, "H", 0); addAtom(aldehyde, "H", 0); // 12 = CHO H
    addBond(aldehyde, 0, 1, SINGLE_BOND, COVALENT_BOND);
    addBond(aldehyde, 1, 2, SINGLE_BOND, COVALENT_BOND);
    addBond(aldehyde, 2, 3, SINGLE_BOND, COVALENT_BOND);
    addBond(aldehyde, 3, 4, DOUBLE_BOND, COVALENT_BOND);
    addBond(aldehyde, 3,12, SINGLE_BOND, COVALENT_BOND);
    addBond(aldehyde, 0, 5, SINGLE_BOND, COVALENT_BOND);
    addBond(aldehyde, 0, 6, SINGLE_BOND, COVALENT_BOND);
    addBond(aldehyde, 0, 7, SINGLE_BOND, COVALENT_BOND);
    addBond(aldehyde, 1, 8, SINGLE_BOND, COVALENT_BOND);
    addBond(aldehyde, 1, 9, SINGLE_BOND, COVALENT_BOND);
    addBond(aldehyde, 2,10, SINGLE_BOND, COVALENT_BOND);
    addBond(aldehyde, 2,11, SINGLE_BOND, COVALENT_BOND);
    printf("isAldehyde : %d\n", isAldehyde(aldehyde));
    printf("Priority   : "); priority(aldehyde);
    deleteMolecule(aldehyde);

    /* ── Amide ── */
    subbanner("Acetamide (Amide)");
    Molecule* amide = createNewMolecule();
    addAtom(amide, "C", 0); // 0 CH3
    addAtom(amide, "C", 0); // 1 carbonyl C
    addAtom(amide, "O", 0); // 2 =O
    addAtom(amide, "N", 0); // 3 NH2
    addAtom(amide, "H", 0); addAtom(amide, "H", 0); addAtom(amide, "H", 0);
    addAtom(amide, "H", 0); addAtom(amide, "H", 0);
    addBond(amide, 0, 1, SINGLE_BOND, COVALENT_BOND);
    addBond(amide, 1, 2, DOUBLE_BOND, COVALENT_BOND);
    addBond(amide, 1, 3, SINGLE_BOND, COVALENT_BOND);
    addBond(amide, 0, 4, SINGLE_BOND, COVALENT_BOND);
    addBond(amide, 0, 5, SINGLE_BOND, COVALENT_BOND);
    addBond(amide, 0, 6, SINGLE_BOND, COVALENT_BOND);
    addBond(amide, 3, 7, SINGLE_BOND, COVALENT_BOND);
    addBond(amide, 3, 8, SINGLE_BOND, COVALENT_BOND);
    printf("isAmide  : %d\n", isAmide(amide));
    printf("Priority : "); priority(amide);
    deleteMolecule(amide);

    /* ── Ketone + Alcohol (FG priority demo) ── */
    subbanner("Hydroxy-ketone (shows FG priority ordering)");
    Molecule* hk = createNewMolecule();
    addAtom(hk, "C", 0); // 0 CH3
    addAtom(hk, "C", 0); // 1 C=O
    addAtom(hk, "C", 0); // 2 CH2
    addAtom(hk, "O", 0); // 3 =O
    addAtom(hk, "O", 0); // 4 OH
    addAtom(hk, "H", 0); addAtom(hk, "H", 0); addAtom(hk, "H", 0);
    addAtom(hk, "H", 0); addAtom(hk, "H", 0); addAtom(hk, "H", 0);
    addBond(hk, 0, 1, SINGLE_BOND, COVALENT_BOND);
    addBond(hk, 1, 2, SINGLE_BOND, COVALENT_BOND);
    addBond(hk, 1, 3, DOUBLE_BOND, COVALENT_BOND);
    addBond(hk, 2, 4, SINGLE_BOND, COVALENT_BOND);
    addBond(hk, 4,10, SINGLE_BOND, COVALENT_BOND);
    addBond(hk, 0, 5, SINGLE_BOND, COVALENT_BOND);
    addBond(hk, 0, 6, SINGLE_BOND, COVALENT_BOND);
    addBond(hk, 0, 7, SINGLE_BOND, COVALENT_BOND);
    addBond(hk, 2, 8, SINGLE_BOND, COVALENT_BOND);
    addBond(hk, 2, 9, SINGLE_BOND, COVALENT_BOND);
    printf("Priority (both OH and CO present): "); priority(hk);
    deleteMolecule(hk);
}

/* ============================================================
   PART 3 : AROMATIC DETECTION
   ============================================================ */
void part3_aromatic() {
    banner("PART 3 : Aromatic Detection (Huckel Rule)");

    /* ── Benzene — should be aromatic ── */
    subbanner("Benzene → expected: aromatic");
    Molecule* benz = createNewMolecule();
    for (int i = 0; i < 6; i++) addAtom(benz, "C", 0);
    for (int i = 0; i < 6; i++) addAtom(benz, "H", 0);
    addBond(benz, 0, 1, DOUBLE_BOND, COVALENT_BOND);
    addBond(benz, 1, 2, SINGLE_BOND,  COVALENT_BOND);
    addBond(benz, 2, 3, DOUBLE_BOND, COVALENT_BOND);
    addBond(benz, 3, 4, SINGLE_BOND,  COVALENT_BOND);
    addBond(benz, 4, 5, DOUBLE_BOND, COVALENT_BOND);
    addBond(benz, 5, 0, SINGLE_BOND,  COVALENT_BOND);
    for (int i = 0; i < 6; i++)
        addBond(benz, i, i + 6, SINGLE_BOND, COVALENT_BOND);
    printf("checkAromatic(Benzene)   : %s\n", checkAromatic(benz) ? "YES — aromatic" : "NO");
    deleteMolecule(benz);

    /* ── Isobutane — no cycle, not aromatic ── */
    subbanner("Isobutane → expected: NOT aromatic (no cycle)");
    Molecule* isob = createNewMolecule();
    addAtom(isob, "C", 0); addAtom(isob, "C", 0);
    addAtom(isob, "C", 0); addAtom(isob, "C", 0);
    for (int i = 0; i < 10; i++) addAtom(isob, "H", 0);
    addBond(isob, 0, 1, SINGLE_BOND, COVALENT_BOND);
    addBond(isob, 0, 2, SINGLE_BOND, COVALENT_BOND);
    addBond(isob, 0, 3, SINGLE_BOND, COVALENT_BOND);
    printf("checkAromatic(Isobutane) : %s\n", checkAromatic(isob) ? "YES" : "NO — not aromatic");
    deleteMolecule(isob);
}

/* ============================================================
   PART 4 : REACTION REPRESENTATION
   ============================================================ */
void part4_reaction() {
    banner("PART 4 : Reaction Representation");

    subbanner("Combustion of Methane: CH4 + O2 → CO2 + H2O");

    // Reactant 1: CH4
    Molecule* ch4 = createNewMolecule();
    addAtom(ch4, "C", 0);
    addAtom(ch4, "H", 0); addAtom(ch4, "H", 0);
    addAtom(ch4, "H", 0); addAtom(ch4, "H", 0);
    addBond(ch4, 0, 1, SINGLE_BOND, COVALENT_BOND);
    addBond(ch4, 0, 2, SINGLE_BOND, COVALENT_BOND);
    addBond(ch4, 0, 3, SINGLE_BOND, COVALENT_BOND);
    addBond(ch4, 0, 4, SINGLE_BOND, COVALENT_BOND);

    // Reactant 2: O2
    Molecule* o2 = createNewMolecule();
    addAtom(o2, "O", 0); addAtom(o2, "O", 0);
    addBond(o2, 0, 1, DOUBLE_BOND, COVALENT_BOND);

    // Product 1: CO2
    Molecule* co2 = createNewMolecule();
    addAtom(co2, "C", 0); addAtom(co2, "O", 0); addAtom(co2, "O", 0);
    addBond(co2, 0, 1, DOUBLE_BOND, COVALENT_BOND);
    addBond(co2, 0, 2, DOUBLE_BOND, COVALENT_BOND);

    // Product 2: H2O
    Molecule* h2o = createNewMolecule();
    addAtom(h2o, "O", 0); addAtom(h2o, "H", 0); addAtom(h2o, "H", 0);
    addBond(h2o, 0, 1, SINGLE_BOND, COVALENT_BOND);
    addBond(h2o, 0, 2, SINGLE_BOND, COVALENT_BOND);

    Reaction* r = createReaction("Combustion of Methane");
    addReactant(r, ch4);
    addReactant(r, o2);
    addProduct(r, co2);
    addProduct(r, h2o);
    printReaction(r);

    // deleteReaction frees the molecules too
    deleteReaction(r);
}

/* ============================================================
   PART 5 : NEW DSA ADDITIONS
   ============================================================ */
void part5_dsa() {
    banner("PART 5 : New DSA Features (dfsTraversal.h)");

    /* helper: build ethanol */
    Molecule* ethanol = createNewMolecule();
    addAtom(ethanol, "C", 0); addAtom(ethanol, "C", 0);
    addAtom(ethanol, "O", 0);
    addAtom(ethanol, "H", 0); addAtom(ethanol, "H", 0);
    addAtom(ethanol, "H", 0); addAtom(ethanol, "H", 0);
    addAtom(ethanol, "H", 0); addAtom(ethanol, "H", 0);
    addBond(ethanol, 0, 1, SINGLE_BOND, COVALENT_BOND);
    addBond(ethanol, 1, 2, SINGLE_BOND, COVALENT_BOND);
    addBond(ethanol, 2, 8, SINGLE_BOND, COVALENT_BOND);
    addBond(ethanol, 0, 3, SINGLE_BOND, COVALENT_BOND);
    addBond(ethanol, 0, 4, SINGLE_BOND, COVALENT_BOND);
    addBond(ethanol, 0, 5, SINGLE_BOND, COVALENT_BOND);
    addBond(ethanol, 1, 6, SINGLE_BOND, COVALENT_BOND);
    addBond(ethanol, 1, 7, SINGLE_BOND, COVALENT_BOND);

    /* helper: build benzene */
    Molecule* benzene = createNewMolecule();
    for (int i = 0; i < 6; i++) addAtom(benzene, "C", 0);
    for (int i = 0; i < 6; i++) addAtom(benzene, "H", 0);
    addBond(benzene, 0, 1, DOUBLE_BOND, COVALENT_BOND);
    addBond(benzene, 1, 2, SINGLE_BOND,  COVALENT_BOND);
    addBond(benzene, 2, 3, DOUBLE_BOND, COVALENT_BOND);
    addBond(benzene, 3, 4, SINGLE_BOND,  COVALENT_BOND);
    addBond(benzene, 4, 5, DOUBLE_BOND, COVALENT_BOND);
    addBond(benzene, 5, 0, SINGLE_BOND,  COVALENT_BOND);
    for (int i = 0; i < 6; i++)
        addBond(benzene, i, i + 6, SINGLE_BOND, COVALENT_BOND);

    /* helper: build isobutane */
    Molecule* isobutane = createNewMolecule();
    addAtom(isobutane, "C", 0); addAtom(isobutane, "C", 0);
    addAtom(isobutane, "C", 0); addAtom(isobutane, "C", 0);
    for (int i = 0; i < 10; i++) addAtom(isobutane, "H", 0);
    addBond(isobutane, 0, 1, SINGLE_BOND, COVALENT_BOND);
    addBond(isobutane, 0, 2, SINGLE_BOND, COVALENT_BOND);
    addBond(isobutane, 0, 3, SINGLE_BOND, COVALENT_BOND);
    addBond(isobutane, 0, 4,  SINGLE_BOND, COVALENT_BOND);
    addBond(isobutane, 1, 5,  SINGLE_BOND, COVALENT_BOND);
    addBond(isobutane, 1, 6,  SINGLE_BOND, COVALENT_BOND);
    addBond(isobutane, 1, 7,  SINGLE_BOND, COVALENT_BOND);
    addBond(isobutane, 2, 8,  SINGLE_BOND, COVALENT_BOND);
    addBond(isobutane, 2, 9,  SINGLE_BOND, COVALENT_BOND);
    addBond(isobutane, 2, 10, SINGLE_BOND, COVALENT_BOND);
    addBond(isobutane, 3, 11, SINGLE_BOND, COVALENT_BOND);
    addBond(isobutane, 3, 12, SINGLE_BOND, COVALENT_BOND);
    addBond(isobutane, 3, 13, SINGLE_BOND, COVALENT_BOND);

    /* 5A — Molecular formula */
    subbanner("5A: Molecular Formula");
    printMolecularFormula(ethanol);
    printMolecularFormula(benzene);
    printMolecularFormula(isobutane);

    /* 5B — Connectivity check */
    subbanner("5B: Connectivity Check (BFS)");
    checkConnectivity(ethanol);
    checkConnectivity(benzene);
    checkConnectivity(isobutane);

    /* 5C — Atom degrees */
    subbanner("5C: Atom Degrees");
    printf("\nEthanol:\n");   printAtomDegrees(ethanol);
    printf("\nBenzene:\n");   printAtomDegrees(benzene);
    printf("\nIsobutane:\n"); printAtomDegrees(isobutane);

    /* 5D — Iterative DFS with manual stack */
    subbanner("5D: Iterative DFS — Manual Stack (no recursion)");
    printf("\n>> Ethanol DFS:");
    dfsIterative(ethanol, 0);

    printf("\n>> Benzene DFS:");
    dfsIterative(benzene, 0);

    printf("\n>> Isobutane DFS:");
    dfsIterative(isobutane, 0);

    deleteMolecule(ethanol);
    deleteMolecule(benzene);
    deleteMolecule(isobutane);
}

/* ============================================================
   MAIN
   ============================================================ */
int main() {
    printf("\n");
    printf("############################################################\n");
    printf("#                                                          #\n");
    printf("#        MOLECULES PROJECT — COMPLETE TEST SUITE          #\n");
    printf("#        Graph-based Molecular Analysis in C              #\n");
    printf("#                                                          #\n");
    printf("############################################################\n");

    part1_core();
    part2_functional_groups();
    part3_aromatic();
    part4_reaction();
    part5_dsa();

    printf("\n");
    printf("############################################################\n");
    printf("#   All parts done. Full project demonstrated.            #\n");
    printf("############################################################\n\n");

    return 0;
}
