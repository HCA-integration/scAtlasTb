ENVS_DIR=$(dirname $0)
TARGET_DIR=/home/icb/chelseaalexandra.bri/miniconda3/envs

for file in "$ENVS_DIR"/*.yaml; do 
    if [ -f "$file" ]; then
        ENV_NAME=$(basename "$file" .yaml)
        SRC_PATH="/lustre/groups/luckylab/workspace/chelsea.bright/envs/$ENV_NAME"   # real env dir
        DEST_PATH="$TARGET_DIR/$ENV_NAME"

        echo "Linking $SRC_PATH → $DEST_PATH"

        # If symlink or directory/file already exists, remove it
        if [ -L "$DEST_PATH" ] || [ -d "$DEST_PATH" ] || [ -f "$DEST_PATH" ]; then
            echo "  Removing existing $DEST_PATH"
            rm -rf "$DEST_PATH"
        fi

        ln -s "$SRC_PATH" "$DEST_PATH"
    fi 
done
