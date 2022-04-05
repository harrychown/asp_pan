# Analyse repeated genes
awk -F "\"" '{ print $4 }' mock_pangenome.gtf | sort | uniq -c | sort -k1 -h | sed 's/^[[:space:]]*//g' | head -n -1 | sed 's/ /,/g' > frequency.csv
grep -v "1," frequency.csv > multifreq.csv
# Identify the frequency of frequencies
awk -F "," '{ print $1 }' frequency.csv | uniq -c | sed 's/^[[:space:]]*//g' | sed 's/ /,/g' > frequency_frequencies.csv
