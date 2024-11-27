

profile_dir=$1

echo " "
ctr="0"
for i in {1..54}; do

    r=$(( 201 + $i ))
    g=$(( (( $i * 4))+39 ))
    b=$(( 54- $i   ))

    rr=$(( 201 + $i ))
        gg=$(( 0 ))
        bb=$(( 147 + $i*2))

    c=$(echo """printf %"$(( $i + 1 ))"s""")
    x=$($c)
    rs=$(echo "$x" | sed 's/ /\\/g')
    if [[ "$ctr" == "2" ]]; then
            colr   "$rs" "$r, $g, $b" "$rr, $gg, $bb" "2"

    fi;
    ctr=$(( $ctr +1 ))
    if [[ "$ctr" == "3" ]]; then
        ctr="1"
    fi
done;

colr -n "  =================||>" "$DY_WT0" "$DY_WB0" "$DY_ES1" && colr -n "  WARNING  " "$DY_ET0" "black" "$DY_ES0" && colr "<||=================  " "$DY_WT0" "$DY_WB0" "$DY_ES1"
colr "  One+ of the template config files in :               " "$DY_ET2" "$DY_EB1" "$DY_ES2"
colr "$profile_dir/templates" "$DY_IT0" "$DY_IB1" "$DY_IS2"
colr "     is newer than your active config files.           " "$DY_ET2" "$DY_EB1" "$DY_ES2"
colr "  You hace 2 options. 1) remove the active files  :    " "$DY_ET2" "$DY_EB1" "$DY_ES2"
colr "rm $profile_dir/*  " "$DY_IT0" "$DY_IB1" "$DY_IS2"

colr "    so the new template files can be copied there. Or  " "$DY_ET2" "$DY_EB1" "$DY_ES2"
colr "    2) Keep the active local copies.  You may touch    " "$DY_ET2" "$DY_EB1" "$DY_ES2"
colr "    them to avoid this error:                          " "$DY_ET2" "$DY_EB1" "$DY_ES2"
colr "touch $profile_dir/*.yaml " "$DY_IT0" "$DY_IB1" "$DY_IS2"
colr "     will keep the files and avoid this warning/block. " "$DY_ET2" "$DY_EB1" "$DY_ES2"
echo " "
echo " "

