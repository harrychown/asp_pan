PREFIX=C30
stats.sh ./mask_$PREFIX/$PREFIX.contigs.fa.masked format=2 > ./$PREFIX.stats
STAT=$(grep '% main' ./$PREFIX.stats | cut -d$'\t' -f2 | sed 's/%//')
MINCONTIG=0
CUTOFF=40
if [ 1 -eq "$(echo "${STAT} < ${CUTOFF}" | bc)" ]
then
	MINCONTIG=5000
else
	MINCONTIG=1000
fi

