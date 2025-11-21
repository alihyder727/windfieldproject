#!/bin/bash

# Check if a directory argument is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <directory>"
    exit 1
fi

# Set the target directory
DIR="$1"

# Iterate over all .gz files in the directory
for file in "$DIR"/*.gz; do
    # Check if the file exists (to prevent errors if no .gz files are present)
    [ -e "$file" ] || continue

    # Skip files that already end with .fits.gz
    if [[ "$file" == *.fits.gz ]]; then
        echo "Skipping: $file (already correctly named)"
        continue
    fi

    # Rename the file to .fits.gz
    new_file="${file%.gz}.fits.gz"
    mv "$file" "$new_file"
    echo "Renamed: $file -> $new_file"
done

echo "Renaming complete."
