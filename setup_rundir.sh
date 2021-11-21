##setup new species run

#make the directory, copy sample sheet, set up run script

SPCODE=$1
SPECIES=$2

mkdir -p $SPCODE
cd $SPCODE
cp ../compPopGen_ms/SRA/cleaned-metadata/${SPECIES}_run_metadata.csv .
git clone https://github.com/harvardinformatics/shortRead_mapping_variantCalling/
mv shortRead_mapping_variantCalling/* .
mv shortRead_mapping_variantCalling/.* .
rmdir shortRead_mapping_variantCalling

#uncomment the following code to switch to the bugfix or dev branch, note newbranch must match something on origin, or comment to use main
#newbranch="bugfix"
#git checkout -b $newbranch
#git push --set-upstream origin $newbranch
#git pull

perl -p -i -e "s/-J sm/-J sm_${SPCODE}/" run_pipeline.sh
perl -p -i -e "s/-o out/-o ${SPCODE}-%j.out/" run_pipeline.sh
perl -p -i -e "s/-e err/-e ${SPCODE}-%j.err/" run_pipeline.sh

#pick one of the following lines depending on coverage
replace="snakemake --snakefile workflow/Snakefile --profile ./profiles/slurm --config samples=\"${SPECIES}_run_metadata.csv\" minD=2 minP=4"
#replace="snakemake --snakefile workflow/Snakefile --profile ./profiles/slurm --config samples=\"${SPECIES}_run_metadata.csv\""

head -n -1 run_pipeline.sh > ${SPCODE}.sh
echo $replace >> ${SPCODE}.sh
