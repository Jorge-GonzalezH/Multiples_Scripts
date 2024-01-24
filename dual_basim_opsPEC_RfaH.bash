# Set the states for RfaH as variables

stateA=opsPEC_RfaH_active
stateB=opsPEC_RfaH_recruited

# Set the first and last atom representing RfaH in the system
first='26931'
last='28225'

# START MESSAGES
echo ' '
echo '### CREATING A DUAL BASIN MODEL - ALL ATOM - GAUSSIAN INTERACTIONS'
echo ' '

# DIHEDRALS
echo '### DIHEDRALS'
echo ' '

# Extract all dihedrals from the top files of each state
awk '{if(($5==1 && NF>7) || ($5==2 && NF>6)) print}' $stateA.top > $stateA.all.dihed
awk '{if(($5==1 && NF>7) || ($5==2 && NF>6)) print}' $stateB.top > $stateB.all.dihed

# Extract only the RfaH dihedrals from each state
awk -v f=$first -v l=$last '{if($1>=f && $4<=l) print}' $stateA.all.dihed > $stateA.dihed
awk -v f=$first -v l=$last '{if($1>=f && $4<=l) print}' $stateB.all.dihed > $stateB.dihed

# Eliminate the RfaH dihedrals from the list of all dihedrals and save in a different file
awk -v f=$first -v l=$last '{if($1>=f && $4<=l) next; print}' $stateA.all.dihed > $stateA.keep.dihed
awk -v f=$first -v l=$last '{if($1>=f && $4<=l) next; print}' $stateB.all.dihed > $stateB.keep.dihed

# Report the number of lines per file
echo 'Total number of dihedrals for $stateA is:'
wc $stateA.all.dihed | awk '{print $1}'
echo ' '
echo 'Total number in the split dihedral files for $stateA is:'
wc $stateA.dihed $stateA.keep.dihed | awk '{print $1, $4}'
echo ' '

echo 'Total number of dihedrals for $stateB is:'
wc $stateB.all.dihed | awk '{print $1}' 
echo ' '
echo 'Total number in the split dihedral files for $stateB is:'
wc $stateB.dihed $stateB.keep.dihed | awk '{print $1, $4}'
echo ' '

echo '### CREATING DUAL DIHEDRAL FILES'

# Generate type 6 dual dihedrals if they diverge by >20 deg between states
paste $stateA.dihed $stateB.dihed | awk '{if($5==2)printf "%-8s%-8s%-8s%-8s%-8s%16.9e%16.9e\n", $1, $2, $3, $4, $5, $6, $7;else if(($8==1 && (($6-$14>20) || ($14-$6>20)))) printf "%-8s%-8s%-8s%-8s%-8s%16.9e%16.9e%16.9e\n", $1, $2, $3, $4, 6, $7, $6, $14; else if(($8==3 && (($6-$14>60) || ($14-$6>60)))) next; else printf "%-8s%-8s%-8s%-8s%-8s%16.9e%16.9e%2s\n", $1, $2, $3, $4, $5, $6, $7, $8}' > $stateA.dual.dihed
paste $stateB.dihed $stateA.dihed | awk '{if($5==2)printf "%-8s%-8s%-8s%-8s%-8s%16.9e%16.9e\n", $1, $2, $3, $4, $5, $6, $7;else if(($8==1 && (($6-$14>20) || ($14-$6>20)))) printf "%-8s%-8s%-8s%-8s%-8s%16.9e%16.9e%16.9e\n", $1, $2, $3, $4, 6, $7, $6, $14; else if(($8==3 && (($6-$14>60) || ($14-$6>60)))) next; else printf "%-8s%-8s%-8s%-8s%-8s%16.9e%16.9e%2s\n", $1, $2, $3, $4, $5, $6, $7, $8}' > $stateB.dual.dihed

# Report the number of dihedrals, counting the dual dihedrals twice
# due to the multiplicity being hardcoded on GROMACS
echo ' '
echo 'Total number of dihedrals in dual dihedral files for $stateA is:'
awk '{if($5==6) printf "%s\n%s\n", $0, $0; else printf "%s\n", $0}' $stateA.dual.dihed | wc | awk '{print $1}'

echo ' '
echo 'Total number of dihedrals in dual dihedral files for $stateB is:'
awk '{if($5==6) printf "%s\n%s\n", $0, $0; else printf "%s\n", $0}' $stateB.dual.dihed | wc | awk '{print $1}'
echo ' '

# Add dual dihedrals into a single file
mv $stateA.dual.dihed tmp
cat $stateA.keep.dihed tmp > $stateA.dual.dihed
rm tmp
mv $stateB.dual.dihed tmp
cat $stateB.keep.dihed tmp > $stateB.dual.dihed
rm tmp

# CONTACT PAIRS
echo '### NATIVE CONTACT PAIRS'
echo ' '

# Extract the list of contacts for RfaH in each state
awk -v f=$first -v l=$last '{if($4>=f && $4<=l) print}' $stateA.contacts.ShadowOutput > $stateA.contacts
awk -v f=$first -v l=$last '{if($4>=f && $4<=l) print}' $stateB.contacts.ShadowOutput > $stateB.contacts

# Extract the list of contacts that will be unmodified
awk -v f=$first -v l=$last '{if($4>=f && $4<=l) next; print}' $stateA.contacts.ShadowOutput > $stateA.keep.contacts
awk -v f=$first -v l=$last '{if($4>=f && $4<=l) next; print}' $stateB.contacts.ShadowOutput > $stateB.keep.contacts

# Extract the contacts that are unique to each state
grep -Fxvf $stateB.contacts $stateA.contacts > $stateA.unique.contacts
grep -Fxvf $stateA.contacts $stateB.contacts > $stateB.unique.contacts

# Extract the contacts that are common to both states
grep -Fxf $stateB.contacts $stateA.contacts > $stateA.common.contacts
grep -Fxf $stateA.contacts $stateB.contacts > $stateB.common.contacts

# Extract all pairs from the top files of each state
awk '{if($3==6 && NF==7) print}' $stateA.top > $stateA.all.pairs
awk '{if($3==6 && NF==7) print}' $stateB.top > $stateB.all.pairs

# Extract pairs with and within RfaH from the pairs files
awk -v f=$first -v l=$last '{if($2>=f && $2<=l) print}' $stateA.all.pairs > $stateA.pairs
awk -v f=$first -v l=$last '{if($2>=f && $2<=l) print}' $stateB.all.pairs > $stateB.pairs

# Extract pairs that will remain unmodified from both systems
awk -v f=$first -v l=$last '{if($2>=f && $2<=l) next; print}' $stateA.all.pairs > $stateA.keep.pairs
awk -v f=$first -v l=$last '{if($2>=f && $2<=l) next; print}' $stateB.all.pairs > $stateB.keep.pairs

# Extract pairs for the common contacts between states
commonrows=`awk 'END{print NR}' $stateA.common.contacts`
if [ -f *.common.pairs ];
then
	rm *.common.pairs
fi

# We read each line in the common contacts list and find the corresponding
# native contact in the pairs files
for ((n=1; n<=commonrows; n++))
do
i=`awk -v n=$n '{if(NR==n) print $2}' $stateA.common.contacts`
j=`awk -v n=$n '{if(NR==n) print $4}' $stateA.common.contacts`
awk -v i=$i -v j=$j '{if($1==i && $2==j) print}' $stateA.pairs >> $stateA.common.pairs
awk -v i=$i -v j=$j '{if($1==i && $2==j) print}' $stateB.pairs >> $stateB.common.pairs
done

# Extract unique pairs to each state into a single pairs file
grep -Fxvf $stateA.common.pairs $stateA.pairs > $stateA.unique.pairs
grep -Fxvf $stateB.common.pairs $stateB.pairs > $stateB.unique.pairs
cat $stateA.keep.pairs $stateA.unique.pairs $stateB.unique.pairs > $stateA.dual.pairs
cat $stateB.keep.pairs $stateA.unique.pairs $stateB.unique.pairs > $stateB.dual.pairs

# Report the total number of contacts for each state
echo 'Total number of contacts for $stateA:'
wc $stateA.all.pairs | awk '{print $1}'
echo ' '
echo 'Number of contacts for RfaH in $stateA:'
wc $stateA.pairs | awk '{print $1}'
echo ' '
echo 'Number of contacts unique to RfaH in $stateA:'
wc $stateA.unique.pairs | awk '{print $1}'
echo ' '

echo 'Total number of contacts for $stateB:'
wc $stateB.all.pairs | awk '{print $1}'
echo ' '
echo 'Number of contacts for RfaH in $stateB:'
wc $stateB.pairs | awk '{print $1}'
echo ' '
echo 'Number of contacts unique to RfaH in $stateB:'
wc $stateB.unique.pairs | awk '{print $1}'
echo ' '

echo 'Number of contacts common to both states'
wc $stateA.common.pairs | awk '{print $1}'
echo ' '


# Generate dual basin contacts for common contacts with distance differences >10%
paste $stateA.common.pairs $stateB.common.pairs | awk '{if(($5/$12)>1.1 || ($12/$5)>1.1) printf "%-8s%-8s%-9s%-16.9e%-16.9e%-16.9e%-16.9e%-16.9e%-16.9e\n", $1, $2, 7, $4, $5, $6, $12, $13, $7; else printf "%-8s%-8s%-9s%-16.9e%-16.9e%-16.9e%-16.9e\n", $1, $2, $3, $4, $5, $6, $7}' >> $stateA.dual.pairs
paste $stateB.common.pairs $stateA.common.pairs | awk '{if(($5/$12)>1.1 || ($12/$5)>1.1) printf "%-8s%-8s%-9s%-16.9e%-16.9e%-16.9e%-16.9e%-16.9e%-16.9e\n", $1, $2, 7, $11, $5, $6, $12, $13, $14; else printf "%-8s%-8s%-9s%-16.9e%-16.9e%-16.9e%-16.9e\n", $8, $9, $10, $11, $12, $13, $14}' >> $stateB.dual.pairs

# Report the dual basin contacts in the final list of pairs
echo 'Total number of dual basin native contacts for $stateA:'
awk -v f=$first -v l=$last '{if($2>=f && $2<=l) print}' $stateA.dual.pairs | wc | awk '{print $1}'

echo ' '
echo 'Total number of dual basin native contacts for $stateB:'
awk -v f=$first -v l=$last '{if($2>=f && $2<=l) print}' $stateB.dual.pairs | wc | awk '{print $1}'
echo ' '

# EXCLUSIONS

echo '### EXCLUSIONS'
echo ' '

# Generate the exclusions lists
awk '{printf "%-8s%-9s\n", $1, $2}' $stateA.dual.pairs > $stateA.dual.exclusions
awk '{printf "%-8s%-9s\n", $1, $2}' $stateB.dual.pairs > $stateB.dual.exclusions

# Report the exclusions for the final topology file
echo 'Total number of native contacts in the exclusion lists for $stateA:'
wc $stateA.dual.exclusions | awk '{print $1}'
echo ' '

echo 'Total number of native contacts in the exclusion lists for $stateB:'
wc $stateB.dual.exclusions | awk '{print $1}'
echo ' '

# Remove all intermediate files
# This can be avoided if required

rm $stateA.contacts $stateB.contacts $stateA.dihed $stateB.dihed *.all.* *.common.* *.keep.* *.unique.*

# GENERATE FINAL TOPOLOGY FILE

echo '### GENERATING DUAL BASIN TOPOLOGY FILE'
echo ' '

# Set several strings that match the structure of the original topology files
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

# Creating the dual topology file for state A
# Copy all text preceding the dihedrals header and add the dual dihedrals
awk -v NRdihedrals=$NRdihedralsA 'NR<=(NRdihedrals+1)' $stateA.top > $stateA.dual.top
paste $stateA.dual.dihed >> $stateA.dual.top

# Copy the pairs header and then add the dual pairs
awk -v NRpairs=$NRpairsA 'NR>=(NRpairs-1) && NR<=(NRpairs+1)' $stateA.top >> $stateA.dual.top
paste $stateA.dual.pairs >> $stateA.dual.top

# Copy the exclusions header and then add the exclusions
awk -v NRexclusions=$NRexclusionsA ' NR>=(NRexclusions-1) && NR<=(NRexclusions+1)' $stateA.top >> $stateA.dual.top
paste $stateA.dual.exclusions >> $stateA.dual.top

# Add the final lines of the original topology file
awk -v NRsystem=$NRsystemA 'NR>=(NRsystem-1)' $stateA.top >> $stateA.dual.top

# Creating the dual topology file for state B
# Copy all text preceding the dihedrals header and add the dual dihedrals
awk -v NRdihedrals=$NRdihedralsB 'NR<=(NRdihedrals+1)' $stateB.top > $stateB.dual.top
paste $stateB.dual.dihed >> $stateB.dual.top

# Copy the pairs header and then add the dual pairs
awk -v NRpairs=$NRpairsB 'NR>=(NRpairs-1) && NR<=(NRpairs+1)' $stateB.top >> $stateB.dual.top
paste $stateB.dual.pairs >> $stateB.dual.top

# Copy the exclusions header and then add the exclusions
awk -v NRexclusions=$NRexclusionsB ' NR>=(NRexclusions-1) && NR<=(NRexclusions+1)' $stateB.top >> $stateB.dual.top
paste $stateB.dual.exclusions >> $stateB.dual.top

# Add the final lines of the original topology file
awk -v NRsystem=$NRsystemB 'NR>=(NRsystem-1)' $stateB.top >> $stateB.dual.top

# Report that the files are done
echo ' '
echo 'Dual basin topology files done!'
