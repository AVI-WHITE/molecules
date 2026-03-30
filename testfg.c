#include<stdio.h>
#include "FGdetection.h"
#include "FGpriority.h"



int main (){

// ALOCOHOL
// Molecule* mol = createNewMolecule();
// // Atoms
// addAtom(mol, "C"); // 0
// addAtom(mol, "C"); // 1
// addAtom(mol, "O"); // 2

// addAtom(mol, "H"); // 3
// addAtom(mol, "H"); // 4
// addAtom(mol, "H"); // 5

// addAtom(mol, "H"); // 6
// addAtom(mol, "H"); // 7

// addAtom(mol, "H"); // 8

// // Bonds
// addBond(mol, 0, 1, SINGLE_BOND, COVALENT_BOND); // C-C
// addBond(mol, 1, 2, SINGLE_BOND, COVALENT_BOND); // C-O

// addBond(mol, 0, 3, SINGLE_BOND, COVALENT_BOND);
// addBond(mol, 0, 4, SINGLE_BOND, COVALENT_BOND);
// addBond(mol, 0, 5, SINGLE_BOND, COVALENT_BOND);

// addBond(mol, 1, 6, SINGLE_BOND, COVALENT_BOND);
// addBond(mol, 1, 7, SINGLE_BOND, COVALENT_BOND);

// addBond(mol, 2, 8, SINGLE_BOND, COVALENT_BOND); // O-H




//CARBOXYLIC ACID

// Molecule* mol = createNewMolecule();

// // Atoms
// addAtom(mol, "C"); // 0 (CH3)
// addAtom(mol, "C"); // 1 (COOH carbon)
// addAtom(mol, "O"); // 2 (=O)
// addAtom(mol, "O"); // 3 (OH)

// addAtom(mol, "H"); // 4
// addAtom(mol, "H"); // 5
// addAtom(mol, "H"); // 6
// addAtom(mol, "H"); // 7

// // Bonds
// addBond(mol, 0, 1, SINGLE_BOND, COVALENT_BOND); // C-C

// addBond(mol, 1, 2, DOUBLE_BOND, COVALENT_BOND); // C=O
// addBond(mol, 1, 3, SINGLE_BOND, COVALENT_BOND); // C-O
// addBond(mol, 3, 7, SINGLE_BOND, COVALENT_BOND); // O-H

// addBond(mol, 0, 4, SINGLE_BOND, COVALENT_BOND);
// addBond(mol, 0, 5, SINGLE_BOND, COVALENT_BOND);
// addBond(mol, 0, 6, SINGLE_BOND, COVALENT_BOND);



//AMIDE

// Molecule* mol = createNewMolecule();

// // Atoms
// addAtom(mol, "C"); // 0 (CH3 carbon)
// addAtom(mol, "C"); // 1 (carbonyl carbon)
// addAtom(mol, "O"); // 2 (=O)
// addAtom(mol, "N"); // 3 (NH2)

// addAtom(mol, "H"); // 4
// addAtom(mol, "H"); // 5
// addAtom(mol, "H"); // 6

// addAtom(mol, "H"); // 7
// addAtom(mol, "H"); // 8

// // Bonds
// addBond(mol, 0, 1, SINGLE_BOND, COVALENT_BOND); // C-C

// addBond(mol, 1, 2, DOUBLE_BOND, COVALENT_BOND); // C=O
// addBond(mol, 1, 3, SINGLE_BOND, COVALENT_BOND); // C-N

// // CH3 hydrogens
// addBond(mol, 0, 4, SINGLE_BOND, COVALENT_BOND);
// addBond(mol, 0, 5, SINGLE_BOND, COVALENT_BOND);
// addBond(mol, 0, 6, SINGLE_BOND, COVALENT_BOND);

// // NH2 hydrogens
// addBond(mol, 3, 7, SINGLE_BOND, COVALENT_BOND);
// addBond(mol, 3, 8, SINGLE_BOND, COVALENT_BOND);
// priority(mol);




// ALDEHDYE

// Molecule* mol = createNewMolecule();

// // Atoms
// addAtom(mol, "C"); // 0
// addAtom(mol, "C"); // 1
// addAtom(mol, "C"); // 2
// addAtom(mol, "C"); // 3 (CHO carbon)

// addAtom(mol, "O"); // 4 (=O)

// addAtom(mol, "H"); // 5
// addAtom(mol, "H"); // 6
// addAtom(mol, "H"); // 7

// addAtom(mol, "H"); // 8
// addAtom(mol, "H"); // 9

// addAtom(mol, "H"); // 10
// addAtom(mol, "H"); // 11

// addAtom(mol, "H"); // 12 (CHO hydrogen)

// // Bonds
// addBond(mol, 0, 1, SINGLE_BOND, COVALENT_BOND);
// addBond(mol, 1, 2, SINGLE_BOND, COVALENT_BOND);
// addBond(mol, 2, 3, SINGLE_BOND, COVALENT_BOND);

// addBond(mol, 3, 4, DOUBLE_BOND, COVALENT_BOND); // C=O
// addBond(mol, 3, 12, SINGLE_BOND, COVALENT_BOND); // C-H (CHO)

// addBond(mol, 0, 5, SINGLE_BOND, COVALENT_BOND);
// addBond(mol, 0, 6, SINGLE_BOND, COVALENT_BOND);
// addBond(mol, 0, 7, SINGLE_BOND, COVALENT_BOND);

// addBond(mol, 1, 8, SINGLE_BOND, COVALENT_BOND);
// addBond(mol, 1, 9, SINGLE_BOND, COVALENT_BOND);

// addBond(mol, 2, 10, SINGLE_BOND, COVALENT_BOND);
// addBond(mol, 2, 11, SINGLE_BOND, COVALENT_BOND);


//FOR FG PRIORITY
Molecule* mol = createNewMolecule();

// Atoms
addAtom(mol, "C"); // 0 (CH3)
addAtom(mol, "C"); // 1 (C=O carbon)
addAtom(mol, "C"); // 2 (CH2)
addAtom(mol, "O"); // 3 (=O)
addAtom(mol, "O"); // 4 (OH)

addAtom(mol, "H"); // 5
addAtom(mol, "H"); // 6
addAtom(mol, "H"); // 7

addAtom(mol, "H"); // 8
addAtom(mol, "H"); // 9

addAtom(mol, "H"); // 10

// Bonds
// main carbon chain
addBond(mol, 0, 1, SINGLE_BOND, COVALENT_BOND);
addBond(mol, 1, 2, SINGLE_BOND, COVALENT_BOND);
// ketone group
addBond(mol, 1, 3, DOUBLE_BOND, COVALENT_BOND); // C=O
// alcohol Group
addBond(mol, 2, 4, SINGLE_BOND, COVALENT_BOND); // C-O
addBond(mol, 4, 10, SINGLE_BOND, COVALENT_BOND); // O-H
// CH3 hydrogens
addBond(mol, 0, 5, SINGLE_BOND, COVALENT_BOND);
addBond(mol, 0, 6, SINGLE_BOND, COVALENT_BOND);
addBond(mol, 0, 7, SINGLE_BOND, COVALENT_BOND);
// CH2 hydrogens
addBond(mol, 2, 8, SINGLE_BOND, COVALENT_BOND);
addBond(mol, 2, 9, SINGLE_BOND, COVALENT_BOND);
priority(mol);





// if(isAlcohol(mol)) printf("Alcohol is detetcted\n");
// else printf("Alcohol is not present\n");

// if(isAldehyde(mol)) printf("Aldehyde is detetcted\n");
// else printf("Aldehyde is not present\n");

// if(isKetone(mol)) printf("Ketone is detetcted\n");
// else printf("Ketone is not present\n");

// if(isCarboxylic(mol)) printf("Carboxylic acid is detetcted\n");
// else printf("Carboxylic acid is not present\n");

// if(isAmide(mol)) printf("Amide is detetcted\n");
// else printf("Amide is not present\n");




//Benzene

// Molecule* benzene = createNewMolecule();
//     addAtom(benzene, "C");
//     addAtom(benzene, "C");
//     addAtom(benzene, "C");
//     addAtom(benzene, "C");
//     addAtom(benzene, "C");
//     addAtom(benzene, "C");

//     addAtom(benzene, "H");
//     addAtom(benzene, "H");
//     addAtom(benzene, "H");
//     addAtom(benzene, "H");
//     addAtom(benzene, "H");
//     addAtom(benzene, "H");

//     addBond(benzene, 0, 1, DOUBLE_BOND, COVALENT_BOND);
//     addBond(benzene, 1, 2, SINGLE_BOND, COVALENT_BOND);
//     addBond(benzene, 2, 3, DOUBLE_BOND, COVALENT_BOND);
//     addBond(benzene, 3, 4, SINGLE_BOND, COVALENT_BOND);
//     addBond(benzene, 4, 5, DOUBLE_BOND, COVALENT_BOND);
//     addBond(benzene, 5, 0, SINGLE_BOND, COVALENT_BOND);

//     addBond(benzene, 0, 6, SINGLE_BOND, COVALENT_BOND);
//     addBond(benzene, 1, 7, SINGLE_BOND, COVALENT_BOND);
//     addBond(benzene, 2, 8, SINGLE_BOND, COVALENT_BOND);
//     addBond(benzene, 3, 9, SINGLE_BOND, COVALENT_BOND);
//     addBond(benzene, 4, 10, SINGLE_BOND, COVALENT_BOND);
//     addBond(benzene, 5, 11, SINGLE_BOND, COVALENT_BOND);

//     if(checkAromatic(benzene) == true) printf("Molecule is aromatic");
//     else printf("Molecule is not aromatic");



    return 0;
}




