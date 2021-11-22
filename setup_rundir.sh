##setup new species run

#make the directory, copy sample sheet, set up run script

####  CHANGE THESE ####
SPCODE=tunAlb #change to species code
SPECIES=Thunnus_albacares #change to species name
branch=bugfix #change to main or dev to run on main or dev
CONDAPATH="/n/holyscratch01/informatics/tsackton/snpArcher_conda/" #set to .snakemake/conda for traditional per-directory location

#pick one depending on coverage parameters
highcov="snakemake --snakefile workflow/Snakefile --profile ./profiles/slurm --use-conda --conda-prefix=\"${CONDAPATH}\" --config samples=\"${SPECIES}_run_metadata.csv\" minD=4 minP=2"
lowcov="snakemake --snakefile workflow/Snakefile --profile ./profiles/slurm --use-conda --conda-prefix=\"${CONDAPATH}\" --config samples=\"${SPECIES}_run_metadata.csv\""

#### CHANGE THIS #####
replace=$lowcov

mkdir -p $SPCODE
cd $SPCODE
cp ../compPopGen_ms/SRA/cleaned-metadata/${SPECIES}_run_metadata.csv .
git clone --branch $branch https://github.com/harvardinformatics/shortRead_mapping_variantCalling/
mv shortRead_mapping_variantCalling/* .
mv shortRead_mapping_variantCalling/.* .
rmdir shortRead_mapping_variantCalling

perl -p -i -e "s/-J sm/-J sm_${SPCODE}/" run_pipeline.sh
perl -p -i -e "s/-o out/-o ${SPCODE}-%j.out/" run_pipeline.sh
perl -p -i -e "s/-e err/-e ${SPCODE}-%j.err/" run_pipeline.sh
perl -p -i -e "s/shared/shared,holy-smokes/" run_pipeline.sh

head -n -1 run_pipeline.sh > ${SPCODE}.sh
echo $replace >> ${SPCODE}.sh
