
while read -r file; do
  localfile="./${file#/}"              # local relative path (strip leading '/')
  fitsfile="${localfile%.gz}.fits.gz"  # target name

  if [[ -f "$localfile" ]]; then
    # If an old fits version exists, remove it first
    if [[ -f "$fitsfile" ]]; then
      echo "⚠️  Removing old: $fitsfile"
      rm -f "$fitsfile"
    fi

    mv "$localfile" "$fitsfile"
    echo "✅ Renamed: $localfile -> $fitsfile"

  else
    echo "❌ Missing: $localfile"
  fi

done < combined_list_gz.txt
