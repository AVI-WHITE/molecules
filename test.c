#include<stdio.h>
#include"molecule.h"

int main(){
    Molecule* benzene = createNewMolecule();
    addAtom(benzene, "C", 0);
    addAtom(benzene, "C", 0);
    addAtom(benzene, "C", 0);
    addAtom(benzene, "C", 0);
    addAtom(benzene, "C", 0);
    addAtom(benzene, "C", 0);

    addAtom(benzene, "H", 0);
    addAtom(benzene, "H", 0);
    addAtom(benzene, "H", 0);
    addAtom(benzene, "H", 0);
    addAtom(benzene, "H", 0);
    addAtom(benzene, "H", 0);

    addBond(benzene, 0, 1, DOUBLE_BOND, COVALENT_BOND);
    addBond(benzene, 1, 2, SINGLE_BOND, COVALENT_BOND);
    addBond(benzene, 2, 3, DOUBLE_BOND, COVALENT_BOND);
    addBond(benzene, 3, 4, SINGLE_BOND, COVALENT_BOND);
    addBond(benzene, 4, 5, DOUBLE_BOND, COVALENT_BOND);
    addBond(benzene, 5, 0, SINGLE_BOND, COVALENT_BOND);

    addBond(benzene, 0, 6, SINGLE_BOND, COVALENT_BOND);
    addBond(benzene, 1, 7, SINGLE_BOND, COVALENT_BOND);
    addBond(benzene, 2, 8, SINGLE_BOND, COVALENT_BOND);
    addBond(benzene, 3, 9, SINGLE_BOND, COVALENT_BOND);
    addBond(benzene, 4, 10, SINGLE_BOND, COVALENT_BOND);
    addBond(benzene, 5, 11, SINGLE_BOND, COVALENT_BOND);

    printMolecule(benzene);
    deleteMolecule(benzene);
    benzene = NULL;

    Molecule* salt = createNewMolecule();

    addAtom(salt, "Na", 1);
    addAtom(salt, "Cl", -1);
    addBond(salt, 0, 1, 1, IONIC_BOND);

    printf("\n\nPrinting structure of common salt NaCl");
    printMolecule(salt);
    deleteMolecule(salt);

    printf("\n\n--- Building Na2SO4 (Sodium Sulfate) ---\n");
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
    printf("\n");


    deleteMolecule(na2so4);

    return 0;
}