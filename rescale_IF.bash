#Set the ratio as the number of bRfaH_CTD contacts over the number of aRfaH_CTD contacts
#In this case, this is 470/261

ratio=1.8

#We finally set a loop for printing 0-100% rescaled interface interactions as in our PLoS article
#We will print them into a new folder called interface_rws

folder=interface_rws
#mkdir $folder

for n in 0.66 0.68 0.69
do
awk -v ratio=$ratio -v n=$n '{printf "%-8s%-8s%-9s%-16.9e%-16.9e%-16.9e%-16.9e\n", $1, $2, $3, $4*ratio*n, $5, $6, $7}' aRfaH_interface.pairs > $folder/aRfaH_interface_rw.$n.pairs

#Then, we put together the files for dualtop and save them into the same folder

cat dualtop_ftype7aRfaH.pairs $folder/aRfaH_interface_rw.$n.pairs > $folder/dualtop_ftype7aRfaH.$n.pairs
cat dualtop_ftype7bRfaH.pairs $folder/aRfaH_interface_rw.$n.pairs > $folder/dualtop_ftype7bRfaH.$n.pairs

#Finally, we make the exclusions lists and save them into the same folder

awk '{print $1, $2}' $folder/dualtop_ftype7aRfaH.$n.pairs > $folder/dualtop_ftype7aRfaH.exclusions
awk '{print $1, $2}' $folder/dualtop_ftype7bRfaH.$n.pairs > $folder/dualtop_ftype7bRfaH.exclusions

#Done
done
