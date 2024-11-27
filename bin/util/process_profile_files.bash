# Verify profile directory exists


profile_dir="$DAY_ROOT/config/day_profiles/$DAY_PROFILE"
if [[ ! -d "$profile_dir" ]]; then
    echo "ERROR: Profile directory '$profile_dir' does not exist."
    echo "Have you run 'day-activate'?"
    exit 3
fi

# Create conda prefix directory if it doesn't exist
#conda_prefix_dir="/fsx/resources/environments/conda/$USER/$(hostname)"
#if [[ ! -d "$conda_prefix_dir" ]]; then
#    echo "Creating Conda prefix directory: $conda_prefix_dir"
#    mkdir -p "$conda_prefix_dir" || {
#        echo "ERROR: Failed to create Conda prefix directory."
#        exit 1
#    }
#    chmod a+wrx "/fsx/resources/environments/conda/$USER"
#    chmod a+wrx "/fsx/resources/environments/container/$USER"
#    chown -R "$USER:$USER" "/fsx/resources/environments/conda/$USER"
#    chown -R "$USER:$USER" "/fsx/resources/environments/container/$USER"
#fi

# Initialize profile configuration files
profile_files=("$profile_dir/rule_config_lowcov.yaml" "$profile_dir/cluster.yaml" "$profile_dir/rule_config.yaml")
templates=("$profile_dir/templates/rule_config_lowcov.yaml" "$profile_dir/templates/cluster.yaml" "$profile_dir/templates/rule_config.yaml")

for i in "${!profile_files[@]}"; do
    if [[ ! -f "${profile_files[i]}" || "${templates[i]}" -nt "${profile_files[i]}" ]]; then
        echo "Copying template files to active config files in $profile_dir"
        cp "$profile_dir/templates/"*yaml "$profile_dir/" || {
            echo "ERROR: Failed to copy template files."
            exit 48
        }
        break
    fi
done

# Replace placeholders in config files
re_safe_mgr=$(echo "$DAY_ROOT" | sed "s/\//\\\\\//g")
sed -i "s/.DAY_ROOT\/c/$re_safe_mgr\/c/g" "$profile_dir/"*yaml
bn=$(basename "$profile_dir" | sed 's/-/\\-/g')
sed -i "s/BASENAME_REGSUB/$bn/g" "$profile_dir/"*yaml
u="$USER"
sed -i "s/USER_REGSUB/$u/g" "$profile_dir/"*yaml
hn="$HOSTNAME"
sed -i "s/HOSTNAME/$hn/g" "$profile_dir/"*yaml

# Unlock snakemake on exit
trap 'snakemake --unlock --profile "$profile_dir" &>> ./unlock_fails.log &' EXIT

# Display help if requested
if [[ "$1" =~ ^(-h|--help|help)$ || -z "$@" ]]; then
    echo "For help, please run: dy-r help"
    exit 2
fi