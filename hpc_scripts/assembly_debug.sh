conda activate bbtools
stats.sh ./mask_$PREFIX/$PREFIX.contigs.fa.masked format=2 > ./$PREFIX.stats
STAT=$(grep '% main' ./$PREFIX.stats | cut -d$'\t' -f2 | sed 's/%//')
MINCONTIG=0
if [ "$STAT" lt 40 ] 
then
	MINCONTIG=5000
else
	MINCONTIG=1000
fi
echo $STAT
echo $MINCONTIG
