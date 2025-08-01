#!/bin/bash
#SBATCH --job-name=custom_ref_cellranger_count
#SBATCH --output=custom_ref%j.out
#SBATCH --error=custom_ref%j.err
#SBATCH  -p cpu-nodes
#SBATCH  -N 1
#SBATCH  -n 16
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your-jphiri@mlw.mw
#SBATCH --account=laboratory

module load Miniconda3/23.9.0-0
module load cellranger/8.0.1

set -e # exit on error

# CONFIGURATION
WORKDIR="$PWD/custom_cellranger_reference"
REF_NAME="custom_human_hiv_covid"
FASTQ_PATH="fastq" #Folder with fastq files
THREADS=8
MEM=64

# Create working directory
mkdir -p "$WORKDIR"
cd "$WORKDIR"

# === STEP 1: Download Human Reference (GRCh38) ===
echo "==> Step 1: Downloading human reference genome (GRCh38)..."
wget -q -O GRCh38.fa.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.genome.fa.gz
wget -q -O GRCh38.gtf.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
gunzip -f GRCh38.fa.gz
gunzip -f GRCh38.gtf.gz

# === STEP 2: Download HIV Genome (AF033819.3) ===
echo="==> Downloading HIV genome..."
curl -O ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz
tar -xzf edirect.tar.gz
export PATH=$PATH:$HOME/edirect
# Download HIV genome (AF033819.3)
efetch -db nuccore -id AF033819.3 -format fasta > HIV.fa
# Rename the FASTA Header (required for CellRanger compatibility)
sed -i '1s/.*/>HIV/' HIV.fa

# Create a custom HIV.gtf file
cat <<EOF > HIV.gtf
HIV    manual    gene        790     2292    .    +    .    gene_id "gag"; gene_name "gag";
HIV    manual    transcript  790     2292    .    +    .    gene_id "gag"; transcript_id "gag.1";
HIV    manual    exon        790     2292    .    +    .    gene_id "gag"; transcript_id "gag.1"; exon_number "1"; exon_id "gag.1.exon1";

HIV    manual    gene        2085    5096    .    +    .    gene_id "pol"; gene_name "pol";
HIV    manual    transcript  2085    5096    .    +    .    gene_id "pol"; transcript_id "pol.1";
HIV    manual    exon        2085    5096    .    +    .    gene_id "pol"; transcript_id "pol.1"; exon_number "1"; exon_id "pol.1.exon1";

HIV    manual    gene        6225    8795    .    +    .    gene_id "env"; gene_name "env";
HIV    manual    transcript  6225    8795    .    +    .    gene_id "env"; transcript_id "env.1";
HIV    manual    exon        6225    8795    .    +    .    gene_id "env"; transcript_id "env.1"; exon_number "1"; exon_id "env.1.exon1";

HIV    manual    gene        5831    6045    .    +    .    gene_id "tat"; gene_name "tat";
HIV    manual    transcript  5831    6045    .    +    .    gene_id "tat"; transcript_id "tat.1";
HIV    manual    exon        5831    6045    .    +    .    gene_id "tat"; transcript_id "tat.1"; exon_number "1"; exon_id "tat.1.exon1";
EOF

# ===STEP 3: Download SARS-CoV-2 Genome (NC_045512.2)
echo "==> Downloading SARS-CoV-2 genome..."
efetch -db nuccore -id NC_045512.2 -format fasta > SARS.fa
# Rename the  SARS GENOME
sed -i '1s/.*/>SARS-CoV-2/' SARS.fa

# Create a custom SARS-Cov-2 .gtf file
cat <<EOF > SARS.gtf
SARS-CoV-2    manual    gene        266     21555    .    +    .    gene_id "ORF1ab"; gene_name "ORF1ab";
SARS-CoV-2    manual    transcript  266     21555    .    +    .    gene_id "ORF1ab"; transcript_id "ORF1ab.1";
SARS-CoV-2    manual    exon        266     21555    .    +    .    gene_id "ORF1ab"; transcript_id "ORF1ab.1"; exon_number "1"; exon_id "ORF1ab.1.exon1";

SARS-CoV-2    manual    gene        21563   25384    .    +    .    gene_id "S"; gene_name "S";
SARS-CoV-2    manual    transcript  21563   25384    .    +    .    gene_id "S"; transcript_id "S.1";
SARS-CoV-2    manual    exon        21563   25384    .    +    .    gene_id "S"; transcript_id "S.1"; exon_number "1"; exon_id "S.1.exon1";

SARS-CoV-2    manual    gene        26245   26472    .    +    .    gene_id "E"; gene_name "E";
SARS-CoV-2    manual    transcript  26245   26472    .    +    .    gene_id "E"; transcript_id "E.1";
SARS-CoV-2    manual    exon        26245   26472    .    +    .    gene_id "E"; transcript_id "E.1"; exon_number "1"; exon_id "E.1.exon1";

SARS-CoV-2    manual    gene        26523   27191    .    +    .    gene_id "M"; gene_name "M";
SARS-CoV-2    manual    transcript  26523   27191    .    +    .    gene_id "M"; transcript_id "M.1";
SARS-CoV-2    manual    exon        26523   27191    .    +    .    gene_id "M"; transcript_id "M.1"; exon_number "1"; exon_id "M.1.exon1";
EOF

# === STEP 4: Concatenate all reference files ===
echo "==> Concatenating FASTA and GTFs..."
cat GRCh38.fa HIV.fa SARS.fa > ${REF_NAME}.fa
cat GRCh38.gtf HIV.gtf SARS.gtf > ${REF_NAME}.gft

# === STEP 5: Build Custome Reference ===
echo "==> Building Cell ranger reference with cellranger mkref..."
cellranger mkref --genome="${REF_NAME}" \
		 --fasta="${REF_NAME}.fa" \
		 --genes="${REF_NAME}.gft" \
		 --nthreads="${THREADS}"

# === STEP6: Run cellranger count on All Samples ===
echo "==> Running cellranger count on all samples in ${FASTQ_PATH}/"
SAMPLES=$(ls ${FASTQ_PATH}/*_R1_001.fastq.gz | sed -E 's/.*\/([^_]+)_.*/\1/' | sort | uniq)

for SAMPLE_ID in $SAMPLES; do
  echo "==> Running cellranger count for: $SAMPLE_ID"
  cellranger count --id="${SAMPLE_ID}" \
                   --transcriptome="${WORKDIR}/${REF_NAME}" \
                   --fastqs="${WORKDIR}/${FASTQ_PATH}" \
                   --sample="${SAMPLE_ID}" \
                   --localcores=${THREADS} \
                   --localmem=${MEM}
done

echo "All done! Output saved in directories per sample under ${WORKDIR}/"
