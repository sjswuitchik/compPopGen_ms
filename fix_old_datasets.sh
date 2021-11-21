#check if new ref genome and old ref genome are the same

SPECIES="Athene_cunicularia"
OLDDIR="/n/holylfs05/LABS/informatics/Lab/holylfs/ashultz/CompPopGen/SPECIES_DATASETS/Acunicularia"

REFGENOME=$(cut -f3,3 -d, ${SPECIES}_run_metadata.csv | sed '1d' | sort | uniq)

#check ref genomes -- assume there an accession assembly stats file like GCF_000247815.1_FicAlb1.5_assembly_stats.txt

OLDGENOME=$(ls $OLDDIR/genome/${REFGENOME}*_assembly_stats.txt)
if [ ${#OLDGENOME[@]} != 1 ]; then
  exit 1
fi

SAMPLES=$(cut -f1,1 -d, ${SPECIES}_run_metadata.csv | sed '1d')

mkdir -p data/genome
mkdir -p data/fastq
for ref in $REFGENOME
do
  mkdir -p results/${SPECIES}/${REFGENOME}/01_mappedReads
  mkdir -p results/${SPECIES}/${REFGENOME}/02_bamSumstats
  mkdir -p results/${SPECIES}/${REFGENOME}/03_gvcfs
  mkdir -p results/${SPECIES}/${REFGENOME}/04_genomicsDB
  mkdir -p results/${SPECIES}/${REFGENOME}/05_vcfs
  mkdir -p results/${SPECIES}/${REFGENOME}/intervalFiles
done

#get interval lists -- these will be under genome/*_interval_lists/
#in order to keep snakemake from recreating the lists we also need a 'fake' bedfile that will just be a list of the chromosomes in the intervals
#need to be moved to results/${SPECIES}/${REFGENOME}/intervalFiles list[n].list
#also will get index of interval number

OLDINT=$(ls $OLDDIR/genome/*_interval_lists/*.interval_list)
INTNUM=0
for interval in $OLDINT;
do
   INTIDTEMP=${interval%.interval_list}
   INTID=${INTIDTEMP##*_}
   cp --preserve=timestamps $interval results/${SPECIES}/${REFGENOME}/intervalFiles/list${INTID}.list
   INTNUM=$((INTNUM+1))
done
cat results/${SPECIES}/${REFGENOME}/intervalFiles/list*.list > results/${SPECIES}/${REFGENOME}/intervalFiles/${REFGENOME}_intervals_fb.bed
#set timestamp on bed file to the timestamp of the first interval file
touch -r results/${SPECIES}/${REFGENOME}/intervalFiles/list1.list results/${SPECIES}/${REFGENOME}/intervalFiles/${REFGENOME}_intervals_fb.bed
echo "Found $INTNUM intervals"

#get bams and gvcfs

for samp in $SAMPLES
do
  cp -v --preserve=timestamps $OLDDIR/dedup/$samp.dedup.sorted.bam results/${SPECIES}/${REFGENOME}/01_mappedReads/${samp}_final.bam
  cp -v --preserve=timestamps $OLDDIR/dedup/$samp.dedup.sorted.bai results/${SPECIES}/${REFGENOME}/01_mappedReads/${samp}_final.bai
  cp -v --preserve=timestamps $OLDDIR/stats/$samp.dedup.metrics.txt results/${SPECIES}/${REFGENOME}/02_bamSumstats/${samp}_dedupMetrics.txt
  cp -v --preserve=timestamps $OLDDIR/stats/$samp.alignment_metrics.txt results/${SPECIES}/${REFGENOME}/02_bamSumstats/${samp}_AlnSumMets.txt
  #intentionally skipping validate.txt as that is fast to regenerate and worth running again
  #also skipping coverage as these seem to be different formats
  mkdir -p results/${SPECIES}/${REFGENOME}/03_gvcfs/$samp
  for i in $(eval echo "{1..$INTNUM}"); do
    cp -v --preserve=timestamps $OLDDIR/gvcf/${samp}.${i}.g.vcf.gz results/${SPECIES}/${REFGENOME}/03_gvcfs/$samp/L${i}.raw.g.vcf.gz
    cp -v --preserve=timestamps $OLDDIR/gvcf/${samp}.${i}.g.vcf.gz.tbi results/${SPECIES}/${REFGENOME}/03_gvcfs/$samp/L${i}.raw.g.vcf.gz.tbi
    touch -r results/${SPECIES}/${REFGENOME}/03_gvcfs/$samp/L${i}.raw.g.vcf.gz.tbi results/${SPECIES}/${REFGENOME}/03_gvcfs/$samp/L${i}.done
  done
done
