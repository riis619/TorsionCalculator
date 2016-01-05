# TorsionCalculator
Used in UCSF Chimera command line.

Input:
mol--molecule given by Protein Data Base structure file
resChainIndicies-select protein backbone chain numbers that user wishes to conform
xyzPoint-point in chimera xyz plane you wish to conform the molecule toward

Function: MorphTo
Corforms protein backbone molecules along acceptable torsion of dihedral angles such that last molecule in resChainIndicies is closest to xyzPoint.
