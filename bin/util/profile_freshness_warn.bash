

profile_dir=$DAY_PROFILE_DIR

echo "a"

if [[ ! -f "$DAY_PROFILE_DIR/rule_config.yaml" && ! -f "$DAY_PROFILE_DIR/config.yaml"  && ! -f "$DAY_PROFILE_DIR/profile_env.sh" ]]; then
    echo "b"    

    colr -n " > >> >>> " "$DY_IB0" "$DY_IT0" "$DY_IS0" 1>&2
    echo "c"
    colr "COPYING TEMPLATE FILES TO ACTIVE CONFIG FILES IN $profile_dir - none detected" "$DY_IT0" "$DY_IB0" "$DY_IS1" 1>&2
    echo "d"
    profile_files=("$profile_dir/rule_config_lowcov.yaml" "$profile_dir/cluster.yaml" "$profile_dir/rule_config.yaml")
    templates=("$profile_dir/templates/rule_config_lowcov.yaml" "$profile_dir/templates/cluster.yaml" "$profile_dir/templates/rule_config.yaml")

    for i in "${!profile_files[@]}"; do
        if [[ ! -f "${profile_files[i]}" || "${templates[i]}" -nt "${profile_files[i]}" ]]; then
            echo "Copying template files to active config files in $profile_dir"
            cp "$profile_dir/templates/"*yaml "$profile_dir/" || {
                echo "ERROR: Failed to copy template files."
                exit 48
            }
            cp "$profile_dir/templates/"*bash "$profile_dir/" || {
                echo "ERROR: Failed to copy template files."
                exit 48
            }

            break
        fi
    done

    echo "e"
    # Replace placeholders in config files
    re_safe_mgr=$(echo "$DAY_ROOT" | sed "s/\//\\\\\//g")
    sed -i "s/.DAY_ROOT\/c/$re_safe_mgr\/c/g" "$profile_dir/"*yaml
    bn=$(basename "$profile_dir" | sed 's/-/\\-/g')
    sed -i "s/BASENAME_REGSUB/$bn/g" "$profile_dir/"*yaml
    u="$USER"
    sed -i "s/USER_REGSUB/$u/g" "$profile_dir/"*yaml
    hn="$HOSTNAME"
    sed -i "s/HOSTNAME/$hn/g" "$profile_dir/"*yaml

elif [[  "config/day_profiles/$DAY_PROFILE/templates/rule_config_lowcov.yaml" -nt  "config/day_profiles/$DAY_PROFILE/rule_config_lowcov.yaml" || "config/day_profiles/$DAY_PROFILE/templates/cluster.yaml" -nt "config/day_profiles/$DAY_PROFILE/cluster.yaml" || "config/day_profiles/$DAY_PROFILE/templates/config.yaml" -nt "config/day_profiles/$DAY_PROFILE/config.yaml" || "config/day_profiles/$DAY_PROFILE/templates/rule_config.yaml" -nt "config/day_profiles/$DAY_PROFILE/rule_config.yaml" || "config/day_profiles/$DAY_PROFILE/templates/profile.info" -nt "config/day_profiles/$DAY_PROFILE/config.yaml" ]]; then


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
    return 9


else
    colr "Your config files in $profile_dir are newer than the templates. clear 2 go." "$DY_WT1" "$DY_WB2" "$DY_WS1 " 1>2
    sleep 0.1
fi