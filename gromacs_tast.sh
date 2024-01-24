#!/bin/bash

# Set some environment variables 
INPUTS_FILES=`pwd` 
# /user/j/jgonzalez/Gromacs/ADP/MD_Alaninedipep/simulation/long
echo "Work directory set to $INPUTS_FILES"
MDP=$INPUTS_FILES/MDP
echo ".mdp files are stored in $MDP"

# Change to the location of your GROMACS installation
GMX=/opt/software/gromacs/5.0.7/bin
GMPI=/opt/software/gromacs/5.0.7-mpi/bin

#start the work to prepare input simulation prior to simulations steps
echo "--- converting pdb file ---"
$GMX/pdb2gmx_mpi -f ala2.pdb -water tip3p -ff amber03 || exit 1 #also try amber99

echo "--- setting box size ---"
$GMX/editconf_mpi -o box.gro -f conf.gro -bt cubic -d 1.2 || exit 1

echo "--- creating water molecules ---"
$GMX/genbox_mpi -o sol.gro -cp box.gro -cs spc216.gro -p topol.top || exit 1

echo "--- adding NaCl in physiologial concentration ---"
$GMX/grompp_mpi -o iongen.tpr -c sol.gro -f em.mdp || exit 1
echo "13\n" 
$GMX/genion_mpi -o ionized.gro -s iongen.tpr -p topol.top -conc 0.15 || exit 1

echo "Preprocessing stage for MD complete."
