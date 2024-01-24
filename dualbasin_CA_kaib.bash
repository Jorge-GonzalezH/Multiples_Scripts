# Names of the fs and gs files

stateA=2fs_g
stateB=gs2_A_fullsym_g

#### ANGLES ####

# Extract angles from .top files for fsKaiB and gsKaiB

awk '{if($4==1 && NF==6) print}' $stateA.top > $stateA.angles
awk '{if($4==1 && NF==6) print}' $stateB.top > $stateB.angles

# We will extract the angles from fsKaiB and gsKaiB that, when compared to the other state,
# differ in more than 10 degrees

paste $stateA.angles $stateB.angles | awk '{if(($5-$11)>10 || ($5-$11)<-10) printf "%-8s%-8s%-8s%-8s%16.9e%16.9e\n", $1, $2, $3, $4, $5, $6}' > $stateA.diff.angles
paste $stateA.angles $stateB.angles | awk '{if(($5-$11)>10 || ($5-$11)<-10) printf "%-8s%-8s%-8s%-8s%16.9e%16.9e\n", $7, $8, $9, $10, $11, $12}' > $stateB.diff.angles

# We will finally merge these selected angles with the native angles from either fsKaiB or gsKaiB

cat $stateA.angles $stateB.diff.angles > mixed_$stateA.angles
cat $stateB.angles $stateA.diff.angles > mixed_$stateB.angles

#### DIHEDRALS ####

# Extract dihedrals from .top files for fsKaiB and gsKaiB
# We will extract dihedrals for residues 1-50 separate from 51-100

awk '{if($5==1 && NF>7 && ($4<=52 || $4>100 && $4<=152)) print}' $stateA.top > $stateA.NTD.dihed
awk '{if($5==1 && NF>7 && ($4>52 && $4<=100 || $4>152)) print}' $stateA.top > $stateA.CTD.dihed
awk '{if($5==1 && NF>7 && ($4<=52 || $4>100 && $4<=152)) print}' $stateB.top > $stateB.NTD.dihed
awk '{if($5==1 && NF>7 && ($4>52 && $4<=100 || $4>152)) print}' $stateB.top > $stateB.CTD.dihed

# Mixed dihedrals (as in Ramirez-Sarmiento et al, PLos Comput Biol 2015)
# We will generate a mixed dihedral file with either fs or gs NTD

cat $stateA.NTD.dihed $stateA.CTD.dihed $stateB.CTD.dihed > mixed_$stateA.dihed
cat $stateB.NTD.dihed $stateB.CTD.dihed $stateA.CTD.dihed > mixed_$stateB.dihed

# Only in case of using dual basin dihedrals. It requires a modified GROMACS version
# Generate dual basin dihedrals (ftype 6) if angles diverge by more than 20 degrees for CTD
# Then put all together with NTD dihedrals from either fsKaiB or gsKaiB

paste $stateA.CTD.dihed $stateB.CTD.dihed | awk '{if(($8==1 && (($6-$14>20) || ($14-$6>20)))) printf "%-8s%-8s%-8s%-8s%-8s%16.9e%16.9e%16.9e\n", $1, $2, $3, $4, 6, $7, $6, $14; else if(($8==3 && (($6-$14>60) || ($14-$6>60)))) next; else printf "%-8s%-8s%-8s%-8s%-8s%16.9e%16.9e%2s\n", $1, $2, $3, $4, $5, $6, $7, $8}' > dual_ft6_$stateA.CTD.dihed
cat $stateA.NTD.dihed dual_ft6_$stateA.CTD.dihed > dual_ft6_$stateA.dihed

paste $stateB.CTD.dihed $stateA.CTD.dihed | awk '{if(($8==1 && (($6-$14>20) || ($14-$6>20)))) printf "%-8s%-8s%-8s%-8s%-8s%16.9e%16.9e%16.9e\n", $1, $2, $3, $4, 6, $7, $6, $14; else if(($8==3 && (($6-$14>60) || ($14-$6>60)))) next; else printf "%-8s%-8s%-8s%-8s%-8s%16.9e%16.9e%2s\n", $1, $2, $3, $4, $5, $6, $7, $8}' > dual_ft6_$stateB.CTD.dihed
cat $stateB.NTD.dihed dual_ft6_$stateB.CTD.dihed > dual_ft6_$stateB.dihed

#### NATIVE CONTACTS ####

# Extract the pairs list from the .top files for fsKaiB and gsKaiB

awk '{if($3==6 && NF==7) print}' $stateA.top > $stateA.pairs
awk '{if($3==6 && NF==7) print}' $stateB.top > $stateB.pairs

# Then, extract protein-protein contacts and monomer contacts from gsKaiB in different files

awk '{if($1!=$3)print}' $stateB.contacts.CG > $stateB.contacts.inter
awk '{if($1==$3)print}' $stateB.contacts.CG > $stateB.contacts.intra

# Since fsKaiB has less contacts than gsKaiB, we will rescale the energies of fsKaiB
# We will do this by increasing the energy of fsKaiB contacts proportional to
# the ratio against the total number of contacts in gsKaiB

NcontA=`wc $stateA.contacts.CG | awk '{print $1}'`
NcontB=`wc $stateB.contacts.CG | awk '{print $1}'`
ratio=`echo "scale=2; $NcontB / $NcontA" | bc`

# We will create a new pairs file for fsKaiB with epsilon(contacts) scaled by our ratio

awk -v ratio=$ratio '{printf "%-8s%-8s%-9s%-16.9e%-16.9e%-16.9e%-16.9e\n", $1, $2, $3, $4*ratio, $5, $6, $7}' $stateA.pairs > $stateA.scaled.pairs

# Now, identify contacts between identical residue pairs in fsKaiB and gsKaiB monomers

perl -ne 'print if ($seen{$_} .= @ARGV) =~ /10$/' $stateA.contacts.CG $stateB.contacts.intra > common.contacts

# We will now retrieve the information about these common contacts from the pairs files
# First we count how many common contacts there are

rows_common=`awk 'END{print NR}' common.contacts`

# Then we read each line in the common contacts list and find the corresponding
# native contacts in the pairs files

for ((n=1; n<=rows_common; n++))
do
i=`awk -v n=$n '{if(NR==n) print $2}' common.contacts`
j=`awk -v n=$n '{if(NR==n) print $4}' common.contacts`

# Note that we are extracting the contact information from both fsKaiB and gsKaiB
# to compare their distances and because we assigned different epsilon(contact)

awk -v i=$i -v j=$j '{if($1==i && $2==j) print}' $stateA.scaled.pairs >> common.$stateA.pairs
awk -v i=$i -v j=$j '{if($1==i && $2==j) print}' $stateB.pairs >> common.$stateB.pairs
done

# We will use a similar loop to split the monomer and dimer contacts of the gsKaiB pairs file
rows_intra=`awk 'END{print NR}' $stateB.contacts.intra`
for ((n=1; n<=rows_intra; n++))
do
i=`awk -v n=$n '{if(NR==n) print $2}' $stateB.contacts.intra`
j=`awk -v n=$n '{if(NR==n) print $4}' $stateB.contacts.intra`
awk -v i=$i -v j=$j '{if($1==i && $2==j) print}' $stateB.pairs >> $stateB.pairs.intra
done

rows_inter=`awk 'END{print NR}' $stateB.contacts.inter`
for ((n=1; n<=rows_inter; n++))
do
i=`awk -v n=$n '{if(NR==n) print $2}' $stateB.contacts.inter`
j=`awk -v n=$n '{if(NR==n) print $4}' $stateB.contacts.inter`
awk -v i=$i -v j=$j '{if($1==i && $2==j) print}' $stateB.pairs >> $stateB.pairs.inter
done

# We will now start creating a new pairs file with dual-basin native contacts
# First, we will add unique contacts from both fsKaiB and gsKaiB

grep -Fxvf common.$stateA.pairs $stateA.scaled.pairs >> dual_ft6.pairs
grep -Fxvf common.$stateB.pairs $stateB.pairs.intra >> dual_ft6.pairs

# We will also add the dimer contacts
paste $stateB.pairs.inter >> dual_ft6.pairs

# Finally, we will make dual gaussians (ftype 7) for common contacts 
# whose distances in fsKaiB and gsKaiB differ by a factor of 1.2.
# Otherwise, we will keep the single gaussian (ftype 6) from either
# fsKaiB or gsKaiB. This is the reason why we are creating 2 files
factor=1.2

paste common.$stateA.pairs common.$stateB.pairs | awk -v factor=$factor '{if(($5/$12)>factor || ($12/$5)>factor) printf "%-8s%-8s%-9s%-16.9e%-16.9e%-16.9e%-16.9e%-16.9e%-16.9e\n", $1, $2, 7, $4, $5, $6, $12, $13, $7; else printf "%-8s%-8s%-9s%-16.9e%-16.9e%-16.9e%-16.9e\n", $1, $2, $3, $4, $5, $6, $7}' > dual_ft7.$stateA.pairs
cat dual_ft6.pairs dual_ft7.$stateA.pairs > dual_ft6_ft7_$stateA.pairs

paste common.$stateA.pairs common.$stateB.pairs | awk -v factor=$factor '{if(($5/$12)>factor || ($12/$5)>factor) printf "%-8s%-8s%-9s%-16.9e%-16.9e%-16.9e%-16.9e%-16.9e%-16.9e\n", $1, $2, 7, $11, $5, $6, $12, $13, $14; else printf "%-8s%-8s%-9s%-16.9e%-16.9e%-16.9e%-16.9e\n", $8, $9, $10, $11, $12, $13, $14}' > dual_ft7.$stateB.pairs
cat dual_ft6.pairs dual_ft7.$stateB.pairs > dual_ft6_ft7_$stateB.pairs

# To finalize, we create the new exclusions lists

awk '{printf "%-8s%-8s\n", $1, $2}' dual_ft6_ft7_$stateA.pairs > dual_$stateA.exclusions
awk '{printf "%-8s%-8s\n", $1, $2}' dual_ft6_ft7_$stateB.pairs > dual_$stateB.exclusions

# Now we will incorporate all this information into new topology files

NRanglesA=`awk '$2=="angles"{print NR}' $stateA.top`
NRanglesB=`awk '$2=="angles"{print NR}' $stateB.top`

NRdihedralsA=`awk '$2=="dihedrals"{print NR}' $stateA.top`
NRdihedralsB=`awk '$2=="dihedrals"{print NR}' $stateB.top`

NRpairsA=`awk '$2=="pairs"{print NR}' $stateA.top`
NRpairsB=`awk '$2=="pairs"{print NR}' $stateB.top`

NRexclusionsA=`awk '$2=="exclusions"{print NR}' $stateA.top`
NRexclusionsB=`awk '$2=="exclusions"{print NR}' $stateB.top`

NRsystemA=`awk '$2=="system"{print NR}' $stateA.top`
NRsystemB=`awk '$2=="system"{print NR}' $stateB.top`

# Creating a new topology file with bonded information (bonds, angles) from fsKaiB
# Using mixed dihedrals
awk -v NRangles=$NRanglesA 'NR<=(NRangles+1)' $stateA.top >> $stateA.dual_ft1dih.top
paste mixed_$stateA.angles >> $stateA.dual_ft1dih.top
awk -v NRdihedrals=$NRdihedralsA 'NR>=(NRdihedrals-1) && NR<=(NRdihedrals+1)' $stateA.top >> $stateA.dual_ft1dih.top
paste mixed_$stateA.dihed >> $stateA.dual_ft1dih.top
awk -v NRpairs=$NRpairsA 'NR>=(NRpairs-1) && NR<=(NRpairs+1)' $stateA.top >> $stateA.dual_ft1dih.top
paste dual_ft6_ft7_$stateA.pairs >> $stateA.dual_ft1dih.top
awk -v NRexclusions=$NRexclusionsA ' NR>=(NRexclusions-1) && NR<=(NRexclusions+1)' $stateA.top >> $stateA.dual_ft1dih.top
paste dual_$stateA.exclusions >> $stateA.dual_ft1dih.top
awk -v NRsystem=$NRsystemA 'NR>=(NRsystem-1)' $stateA.top >> $stateA.dual_ft1dih.top

# Creating a new topology file with bonded information (bonds, angles) from fsKaiB
# Using dual dihedrals (ftype 6)
awk -v NRangles=$NRanglesA 'NR<=(NRangles+1)' $stateA.top >> $stateA.dual_ft6dih.top
paste mixed_$stateA.angles >> $stateA.dual_ft6dih.top
awk -v NRdihedrals=$NRdihedralsA 'NR>=(NRdihedrals-1) && NR<=(NRdihedrals+1)' $stateA.top >> $stateA.dual_ft6dih.top
paste dual_ft6_$stateA.dihed >> $stateA.dual_ft6dih.top
awk -v NRpairs=$NRpairsA 'NR>=(NRpairs-1) && NR<=(NRpairs+1)' $stateA.top >> $stateA.dual_ft6dih.top
paste dual_ft6_ft7_$stateA.pairs >> $stateA.dual_ft6dih.top
awk -v NRexclusions=$NRexclusionsA ' NR>=(NRexclusions-1) && NR<=(NRexclusions+1)' $stateA.top >> $stateA.dual_ft6dih.top
paste dual_$stateA.exclusions >> $stateA.dual_ft6dih.top
awk -v NRsystem=$NRsystemA 'NR>=(NRsystem-1)' $stateA.top >> $stateA.dual_ft6dih.top

# Creating a new topology file with bonded information (bonds, angles) from gsKaiB
# Using mixed dihedrals
awk -v NRangles=$NRanglesB 'NR<=(NRangles+1)' $stateB.top >> $stateB.dual_ft1dih.top
paste mixed_$stateB.angles >> $stateB.dual_ft1dih.top
awk -v NRdihedrals=$NRdihedralsB 'NR>=(NRdihedrals-1) && NR<=(NRdihedrals+1)' $stateB.top >> $stateB.dual_ft1dih.top
paste mixed_$stateB.dihed >> $stateB.dual_ft1dih.top
awk -v NRpairs=$NRpairsB 'NR>=(NRpairs-1) && NR<=(NRpairs+1)' $stateB.top >> $stateB.dual_ft1dih.top
paste dual_ft6_ft7_$stateB.pairs >> $stateB.dual_ft1dih.top
awk -v NRexclusions=$NRexclusionsB ' NR>=(NRexclusions-1) && NR<=(NRexclusions+1)' $stateB.top >> $stateB.dual_ft1dih.top
paste dual_$stateB.exclusions >> $stateB.dual_ft1dih.top
awk -v NRsystem=$NRsystemB 'NR>=(NRsystem-1)' $stateB.top >> $stateB.dual_ft1dih.top

# Creating a new topology file with bonded information (bonds, angles) from gsKaiB
# Using dual dihedrals (ftype 6)
awk -v NRangles=$NRanglesB 'NR<=(NRangles+1)' $stateB.top >> $stateB.dual_ft6dih.top
paste mixed_$stateB.angles >> $stateB.dual_ft6dih.top
awk -v NRdihedrals=$NRdihedralsB 'NR>=(NRdihedrals-1) && NR<=(NRdihedrals+1)' $stateB.top >> $stateB.dual_ft6dih.top
paste dual_ft6_$stateB.dihed >> $stateB.dual_ft6dih.top
awk -v NRpairs=$NRpairsB 'NR>=(NRpairs-1) && NR<=(NRpairs+1)' $stateB.top >> $stateB.dual_ft6dih.top
paste dual_ft6_ft7_$stateB.pairs >> $stateB.dual_ft6dih.top
awk -v NRexclusions=$NRexclusionsB ' NR>=(NRexclusions-1) && NR<=(NRexclusions+1)' $stateB.top >> $stateB.dual_ft6dih.top
paste dual_$stateB.exclusions >> $stateB.dual_ft6dih.top
awk -v NRsystem=$NRsystemB 'NR>=(NRsystem-1)' $stateB.top >> $stateB.dual_ft6dih.top
