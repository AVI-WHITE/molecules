#include<stdio.h>
#include"molecule.h"
// #include"longestChain.h"
#include"operationsOnMolecules.h"

int main(){

    // ================= BENZENE =================
    printf("\n--- Benzene (C6H6) ---\n");

    Molecule* benzene = createNewMolecule();

    // Add 6 Carbon atoms
    for(int i = 0; i < 6; i++)
        addAtom(benzene, "C", 0);

    // Add 6 Hydrogen atoms
    for(int i = 0; i < 6; i++)
        addAtom(benzene, "H", 0);

    // Carbon ring
    addBond(benzene, 0, 1, DOUBLE_BOND, COVALENT_BOND);
    addBond(benzene, 1, 2, SINGLE_BOND, COVALENT_BOND);
    addBond(benzene, 2, 3, DOUBLE_BOND, COVALENT_BOND);
    addBond(benzene, 3, 4, SINGLE_BOND, COVALENT_BOND);
    addBond(benzene, 4, 5, DOUBLE_BOND, COVALENT_BOND);
    addBond(benzene, 5, 0, SINGLE_BOND, COVALENT_BOND);

    // Hydrogens
    for(int i = 0; i < 6; i++)
        addBond(benzene, i, i+6, SINGLE_BOND, COVALENT_BOND);

    printMolecule(benzene);

    int longest = findLongestCarbonChain(benzene);
    printf("\n\nLongest Carbon Chain Length (Benzene): %d", longest);
    printf("\nIs benzene cyclic : %d", doesMoleculeContainCycle(benzene));

    deleteMolecule(benzene);


    // ================= Na2SO4 =================
    printf("\n\n--- Na2SO4 ---\n");

    Molecule* na2so4 = createNewMolecule();

    addAtom(na2so4, "S", 0);  
    addAtom(na2so4, "O", 0);  
    addAtom(na2so4, "O", 0);  
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

    longest = findLongestCarbonChain(na2so4);
    printf("\n\nLongest Carbon Chain Length (Na2SO4): %d", longest);
    printf("\nMolecular weight : %lf g/mol\n", molecularWeight(na2so4));
    int* oxsNa2SO4 = calculateOxidationStates(na2so4);
    printOxidationStatesArray(na2so4, oxsNa2SO4);

    deleteMolecule(na2so4);
    free(oxsNa2SO4);

    
    
    // ================= Isobutane =================
    Molecule* isobutane = createNewMolecule();
    printf("\n\n--- Isobutane ---\n");

    addAtom(isobutane, "C", 0); // 0 (Central)
    addAtom(isobutane, "C", 0); // 1 (Branch)
    addAtom(isobutane, "C", 0); // 2 (Branch)
    addAtom(isobutane, "C", 0); // 3 (Branch)
    
    // Hydrogens (Indices 4 to 13)
    for(int i = 0; i < 10; i++) {
        addAtom(isobutane, "H", 0); 
    }
    
    // Carbon-Carbon Bonds
    addBond(isobutane, 0, 1, SINGLE_BOND, COVALENT_BOND);
    addBond(isobutane, 0, 2, SINGLE_BOND, COVALENT_BOND);
    addBond(isobutane, 0, 3, SINGLE_BOND, COVALENT_BOND);
    
    // Carbon-Hydrogen Bonds
    addBond(isobutane, 0, 4, SINGLE_BOND, COVALENT_BOND); // Central C gets 1 H
    
    addBond(isobutane, 1, 5, SINGLE_BOND, COVALENT_BOND); // Branch 1 gets 3 H
    addBond(isobutane, 1, 6, SINGLE_BOND, COVALENT_BOND);
    addBond(isobutane, 1, 7, SINGLE_BOND, COVALENT_BOND);
    
    addBond(isobutane, 2, 8, SINGLE_BOND, COVALENT_BOND); // Branch 2 gets 3 H
    addBond(isobutane, 2, 9, SINGLE_BOND, COVALENT_BOND);
    addBond(isobutane, 2, 10, SINGLE_BOND, COVALENT_BOND);
    
    addBond(isobutane, 3, 11, SINGLE_BOND, COVALENT_BOND); // Branch 3 gets 3 H
    addBond(isobutane, 3, 12, SINGLE_BOND, COVALENT_BOND);
    addBond(isobutane, 3, 13, SINGLE_BOND, COVALENT_BOND);

    printf("\nIs isobutane cyclic : %d", doesMoleculeContainCycle(isobutane));
    
    printf("\nLength of longest carbon chain : %d", findLongestCarbonChain(isobutane));
    
    printf("\nMolecular weight : %lf g/mol\n", molecularWeight(isobutane));
    
    int* oxsIsobutane = calculateOxidationStates(isobutane);
    printOxidationStatesArray(isobutane, oxsIsobutane);
    
    free(oxsIsobutane);
    oxsIsobutane = NULL;
    deleteMolecule(isobutane);



    // ================= Naphthalene =================
    Molecule* naphthalene = createNewMolecule();
    printf("\n\n--- Naphthalene ---\n");

    
    // Carbons (Indices 0 to 9)
    for(int i = 0; i < 10; i++) {
        addAtom(naphthalene, "C", 0); 
    }
    
    // Hydrogens (Indices 10 to 17)
    for(int i = 0; i < 8; i++) {
        addAtom(naphthalene, "H", 0); 
    }
    
    // Carbon Skeleton
    addBond(naphthalene, 0, 1, DOUBLE_BOND, COVALENT_BOND);
    addBond(naphthalene, 1, 2, SINGLE_BOND, COVALENT_BOND);
    addBond(naphthalene, 2, 3, DOUBLE_BOND, COVALENT_BOND);
    addBond(naphthalene, 3, 4, SINGLE_BOND, COVALENT_BOND);
    addBond(naphthalene, 4, 9, DOUBLE_BOND, COVALENT_BOND); // Bridge
    addBond(naphthalene, 9, 0, SINGLE_BOND, COVALENT_BOND);
    
    addBond(naphthalene, 4, 5, SINGLE_BOND, COVALENT_BOND);
    addBond(naphthalene, 5, 6, DOUBLE_BOND, COVALENT_BOND);
    addBond(naphthalene, 6, 7, SINGLE_BOND, COVALENT_BOND);
    addBond(naphthalene, 7, 8, DOUBLE_BOND, COVALENT_BOND);
    addBond(naphthalene, 8, 9, SINGLE_BOND, COVALENT_BOND); // Bridge
    
    // Carbon-Hydrogen Bonds (Bridgehead carbons 4 and 9 get NO hydrogens)
    addBond(naphthalene, 0, 10, SINGLE_BOND, COVALENT_BOND);
    addBond(naphthalene, 1, 11, SINGLE_BOND, COVALENT_BOND);
    addBond(naphthalene, 2, 12, SINGLE_BOND, COVALENT_BOND);
    addBond(naphthalene, 3, 13, SINGLE_BOND, COVALENT_BOND);
    
    addBond(naphthalene, 5, 14, SINGLE_BOND, COVALENT_BOND);
    addBond(naphthalene, 6, 15, SINGLE_BOND, COVALENT_BOND);
    addBond(naphthalene, 7, 16, SINGLE_BOND, COVALENT_BOND);
    addBond(naphthalene, 8, 17, SINGLE_BOND, COVALENT_BOND);

    printf("\nIs naphthalene cyclic : %d", doesMoleculeContainCycle(naphthalene));
    
    printf("\nLength of longest carbon chain : %d", findLongestCarbonChain(naphthalene));
    
    printf("\nMolecular weight : %lf g/mol\n", molecularWeight(naphthalene));
    
    int* oxsNaphthalene = calculateOxidationStates(naphthalene);
    printOxidationStatesArray(naphthalene, oxsNaphthalene);
    
    free(oxsNaphthalene);
    oxsIsobutane = NULL;
    deleteMolecule(naphthalene);

    return 0;
}