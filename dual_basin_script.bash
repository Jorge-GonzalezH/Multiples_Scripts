#Separate dihedrals from aRfaH and bRfaH into NTD and non-NTD
awk '{if($4<784) print}' aRfaH.dihed > aRfaH_NTD.dihed
awk '{if($4<784) print}' bRfaH.dihed > bRfaH_NTD.dihed
awk '{if($4>783) print}' aRfaH.dihed > aRfaH_CTD.dihed
awk '{if($4>783) print}' bRfaH.dihed > bRfaH_CTD.dihed

#Generate dihedral with ftype 6 for dual dihedrals if angles diverge by more than 20 degrees for CTD
paste aRfaH_CTD.dihed bRfaH_CTD.dihed | awk '{if($5==2)printf "%-8s%-8s%-8s%-8s%-8s%16.9e%16.9e\n", $1, $2, $3, $4, $5, $6, $7;else if(($8==1 && (($6-$14>20) || ($14-$6>20)))) printf "%-8s%-8s%-8s%-8s%-8s%16.9e%16.9e%16.9e\n", $1, $2, $3, $4, 6, $7, $6, $14; else if(($8==3 && (($6-$14>60) || ($14-$6>60)))) next; else printf "%-8s%-8s%-8s%-8s%-8s%16.9e%16.9e%2s\n", $1, $2, $3, $4, $5, $6, $7, $8}' > dualtop_ftype6_aRfaH_CTD.dihed
paste bRfaH_CTD.dihed aRfaH_CTD.dihed | awk '{if($5==2)printf "%-8s%-8s%-8s%-8s%-8s%16.9e%16.9e\n", $1, $2, $3, $4, $5, $6, $7;else if(($8==1 && (($6-$14>20) || ($14-$6>20)))) printf "%-8s%-8s%-8s%-8s%-8s%16.9e%16.9e%16.9e\n", $1, $2, $3, $4, 6, $7, $6, $14; else if(($8==3 && (($6-$14>60) || ($14-$6>60)))) next; else printf "%-8s%-8s%-8s%-8s%-8s%16.9e%16.9e%2s\n", $1, $2, $3, $4, $5, $6, $7, $8}' > dualtop_ftype6_bRfaH_CTD.dihed

#Lastly, put NTD and CTD dihedrals together
cat aRfaH_NTD.dihed dualtop_ftype6_aRfaH_CTD.dihed > dualtop_ftype6_aRfaH.dihed
cat bRfaH_NTD.dihed dualtop_ftype6_bRfaH_CTD.dihed > dualtop_ftype6_bRfaH.dihed

#Extract common contacts
perl -ne 'print if ($seen{$_} .= @ARGV) =~ /10$/' aRfaH_CTD.contacts bRfaH_CTD.contacts > common_CTD.contacts

#Extract gauss pairs of these common contacts

#First we count how many common contacts there are
rows=`awk 'END{print NR}' common_CTD.contacts`

#Then we extract
for ((n=1; n<=rows; n++))
do
i=`awk -v n=$n '{if(NR==n) print $2}' common_CTD.contacts`
j=`awk -v n=$n '{if(NR==n) print $4}' common_CTD.contacts`

#We extract from both aRfaH and bRfaH to compare distances and due to the different weights
awk -v i=$i -v j=$j '{if($1==i && $2==j) print}' aRfaH_CTD.pairs >> common_aRfaH_CTD.pairs
awk -v i=$i -v j=$j '{if($1==i && $2==j) print}' bRfaH_CTD.pairs >> common_bRfaH_CTD.pairs
done

#We create a new contact file in which we first deposit the NTD contacts
paste NTD.pairs > dualtop.pairs

#We discard the common contacts from aRfaH and bRfaH to add aRfaH- and bRfaH-exclusive contacts
#to our new contact file
grep -Fxvf common_aRfaH_CTD.pairs aRfaH_CTD.pairs >> dualtop.pairs
grep -Fxvf common_bRfaH_CTD.pairs bRfaH_CTD.pairs >> dualtop.pairs

#Finally, we make dual gaussians for those common contacts whose distances
#between aRfaH and bRfaH are different by a factor larger than 1.1
paste common_aRfaH_CTD.pairs common_bRfaH_CTD.pairs | awk '{if(($5/$12)>1.1 || ($12/$5)>1.1) printf "%-8s%-8s%-9s%-16.9e%-16.9e%-16.9e%-16.9e%-16.9e%-16.9e\n", $1, $2, 7, $4, $5, $6, $12, $13, $7; else printf "%-8s%-8s%-9s%-16.9e%-16.9e%-16.9e%-16.9e\n", $1, $2, $3, $4, $5, $6, $7}' > ftype7_aRfaH.pairs
paste common_aRfaH_CTD.pairs common_bRfaH_CTD.pairs | awk '{if(($5/$12)>1.1 || ($12/$5)>1.1) printf "%-8s%-8s%-9s%-16.9e%-16.9e%-16.9e%-16.9e%-16.9e%-16.9e\n", $1, $2, 7, $11, $5, $6, $12, $13, $14; else printf "%-8s%-8s%-9s%-16.9e%-16.9e%-16.9e%-16.9e\n", $8, $9, $10, $11, $12, $13, $14}' > ftype7_bRfaH.pairs
cat dualtop.pairs ftype7_aRfaH.pairs > dualtop_ftype7aRfaH.pairs
cat dualtop.pairs ftype7_bRfaH.pairs > dualtop_ftype7bRfaH.pairs

#Create exclusions list and add to appropriate top
awk '{print $1, $2}' dualtop_ftype7aRfaH.pairs > dualtop_ftype7aRfaH.exclusions
awk '{print $1, $2}' dualtop_ftype7bRfaH.pairs > dualtop_ftype7bRfaH.exclusions

#We will finally have to add the interface contacts but we need to rescale first. We do this separately
