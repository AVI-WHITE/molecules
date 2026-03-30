#include<stdio.h>
#include"molecule.h"
#include"longestChain.h"
#include"reaction.h"


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
    printf("\nLongest Carbon Chain Length (Benzene): %d\n", longest);

    deleteMolecule(benzene);


    // ================= NaCl =================
    printf("\n--- NaCl ---\n");

    Molecule* salt = createNewMolecule();

    addAtom(salt, "Na", 1);
    addAtom(salt, "Cl", -1);
    addBond(salt, 0, 1, 1, IONIC_BOND);

    printMolecule(salt);

    longest = findLongestCarbonChain(salt);
    printf("\nLongest Carbon Chain Length (NaCl): %d\n", longest);

    deleteMolecule(salt);


    // ================= Na2SO4 =================
    printf("\n--- Na2SO4 ---\n");

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
    printf("\nLongest Carbon Chain Length (Na2SO4): %d\n", longest);

    deleteMolecule(na2so4);


    Reaction* r = createReaction("Sample Reaction");

// Step 1: Create molecules
    Molecule* m1 = createNewMolecule();
    addAtom(m1, "C", 0);

    Molecule* m2 = createNewMolecule();
    addAtom(m2, "O", 0);

    

    // Step 2: ADD them to reaction
    addReactant(r, m1);
    addProduct(r, m2);

    // Step 3: Print
    printReaction(r);

    return 0;
}