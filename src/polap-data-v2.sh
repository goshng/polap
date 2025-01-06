#!/usr/bin/bash

# src/polap.sh x-get-sra-info --sra <SRA>

S=(
	SRR14883332
	SRR18079816
	SRR13238610
	ERR338629
	SRR11547303
	SRR1980665
	SRR7153095
	SRR10948618
	SRR8989349
	SRR8989348
	SRR8989347
	SRR8989346
	SRR8989345
	SRR8989344
	SRR21038753
)

for i in "${S[@]}"; do
	# src/polap.sh x-ncbi-fetch-sra --sra "$i"
	# sleep 5
	rm -rf "$i"
done

exit
#!/usr/bin/bash

# src/polap.sh x-get-sra-info --sra <SRA>

S=(SRR7153095
	SRR7161123
	SRR21976089
	SRR21976090
	SRR21976091
	SRR21976092
	SRR14298760
	SRR14298751
	SRR14298746
	SRR14298745

	SRR14883332
	SRR18079816
	SRR13238610
	ERR338629
	SRR11547303
	SRR1980665
	SRR7153095
	SRR10948618
	SRR8989349
	SRR8989348
	SRR8989347
	SRR8989346
	SRR8989345
	SRR8989344
	SRR21038753)

for i in "${S[@]}"; do
	src/polap.sh x-get-sra-info --sra "$i"
	sleep 5
done

exit

# Table 2 of Zhou2023
SRR14298760 35261850229 WGS GENOMIC SINGLE OXFORD_NANOPORE Juncus effusus
SRR14298746 28996069500 WGS GENOMIC PAIRED ILLUMINA Juncus effusus

SRR14298751 37591299287 WGS GENOMIC SINGLE OXFORD_NANOPORE Juncus inflexus
SRR14298745 25023621900 WGS GENOMIC PAIRED ILLUMINA Juncus inflexus

SRR21976090 2302422897 WGS GENOMIC SINGLE OXFORD_NANOPORE Juncus roemerianus
SRR21976092 79461161000 WGS GENOMIC PAIRED ILLUMINA Juncus roemerianus

SRR21976089 561799410 WGS GENOMIC SINGLE OXFORD_NANOPORE Juncus validus
SRR21976091 78356215000 WGS GENOMIC PAIRED ILLUMINA Juncus validus

# Table 1 of Zhou2023: with short-read data as well
SRR7153095 6040851677 WGS GENOMIC SINGLE OXFORD_NANOPORE Eucalyptus pauciflora
SRR7161123 6376665372 WGS GENOMIC PAIRED ILLUMINA Eucalyptus pauciflora

# Table 1 of Zhou2023
SRR14883332 27304352670 WGS GENOMIC SINGLE PACBIO_SMRT Arctostaphylos glauca
SRR18079816 2615947808 WGA GENOMIC SINGLE PACBIO_SMRT Lepidium sativum
SRR13238610 792414508 Targeted-Capture OTHER SINGLE PACBIO_SMRT Chaetoceros muellerii
ERR338629 54492250 WGS GENOMIC SINGLE PACBIO_SMRT Potentilla micrantha
SRR11547303 5987567965 OTHER GENOMIC SINGLE PACBIO_SMRT Durio zibethinus
SRR1980665 296752589 WGS GENOMIC SINGLE PACBIO_SMRT Beta vulgaris subsp. vulgaris
# Eleocharis dulcis - available by asking authors (Ref: NC_047447.1)
SRR7153095 6040851677 WGS GENOMIC SINGLE OXFORD_NANOPORE Eucalyptus pauciflora
SRR10948618 66965150 AMPLICON GENOMIC SINGLE OXFORD_NANOPORE Leucanthemum vulgare
SRR8989349 317851242 Targeted-Capture GENOMIC SINGLE OXFORD_NANOPORE Oryza glaberrima
SRR8989348 525852568 Targeted-Capture GENOMIC SINGLE OXFORD_NANOPORE Cenchrus americanus
SRR8989347 555524185 Targeted-Capture GENOMIC SINGLE OXFORD_NANOPORE Digitaria exilis
SRR8989346 645897952 Targeted-Capture GENOMIC SINGLE OXFORD_NANOPORE Podococcus acaulis
SRR8989345 202292169 Targeted-Capture GENOMIC SINGLE OXFORD_NANOPORE Raphia textilis
SRR8989344 505797880 Targeted-Capture GENOMIC SINGLE OXFORD_NANOPORE Phytelephas aequatorialis
SRR21038753 5928239580 WGS GENOMIC SINGLE OXFORD_NANOPORE Picea glauca

# Eucalyptus pauciflora
src/polap.sh x-ncbi-fetch-sra --sra SRR7153095
src/polap.sh x-ncbi-fetch-sra --sra SRR7161123

#
src/polap.sh x-ncbi-fetch-sra --sra SRR21976089
src/polap.sh x-ncbi-fetch-sra --sra SRR21976090

src/polap.sh x-ncbi-fetch-sra --sra SRR21976091
# J. validus (SRR21976091)
# J. roamerianus (SRR21976092)
src/polap.sh x-ncbi-fetch-sra --sra SRR21976092

#
src/polap.sh x-ncbi-fetch-sra --sra SRR14298760
src/polap.sh x-ncbi-fetch-sra --sra SRR14298751

src/polap.sh x-ncbi-fetch-sra --sra SRR14298746
src/polap.sh x-ncbi-fetch-sra --sra SRR14298745

src/polap.sh x-ncbi-fetch-sra-runinfo --sra "$long_sra"

exit

#Download Long read data or Illumina data
module add sratoolkit/3.0.0

fasterq-dump --split-files SRR14883332
#Download reference genome
#Long read data
esearch -db nucleotide -query "NC_035584.1" | efetch -format fasta >NC_035584.1.fasta

# Replace SRA accession numbers (SRR14883332) and reference number (NC_035584.1) with the following numbers.
# Long read data for ptGAUL validation:

# Arctostaphylos glauca (SRR14883332) (Ref: NC_035584.1/NC_042713.1/NC_047438.1)
SRR14883332

# Lepidium sativum (SRR18079816) (Ref: NC_047178.1)
SRR18079816

# Chaetoceros muellerii (SRR13238610) (Ref: MW004650.1)
SRR13238610

# Potentilla micrantha (ERR338629) (Ref: NC_015206.1)
ERR338629

# Durio zibethinus (SRR11547303) (Ref: MT321069)
SRR11547303
# Beta vulgaris (SRR1980665) (Ref: KR230391.1)
SRR1980665
# Eleocharis dulcis - available by asking authors (Ref: NC_047447.1)

# Eucalyptus pauciflora (SRR7153095) (Ref: MZ670598.1/HM347959.1/NC_014570.1/AY780259.1/ NC_039597.1)
SRR7153095

# Leucanthemum vulgare (SRR10948618) (Ref: NC_047460.1)
SRR10948618

# Oryza glaberrima (SRR8989349) (Ref: NC_024175.1)
SRR8989349

# Cenchrus americanus (SRR8989348) (Ref: NC_024171.1)
SRR8989348

# Digitaria exilis (SRR8989347) (Ref: NC_024176.1)
SRR8989347

# Podococcus acaulis (SRR8989346) (Ref: NC_027276.1)
SRR8989346

# Raphia textilis (SRR8989345) (Ref: NC_020365.1)
SRR8989345

# Phytelephas aequatorialis (SRR8989344) (Ref: NC_029957.1)
SRR8989344

# Picea glauca (SRR21038753) (Ref: NC_021456.1).
SRR21038753

# Illumina data of Juncus for GetOrganelle:

# J. inflexus (SRR14298745)
# J. effusus (SRR13309655 (Lu et al. 2021) and SRR14298746)
# J. roamerianus (SRR21976092)
# J. validus (SRR21976091)

# example of using GetOrganelle for J. effusus
fasterq-dump --split-files SRR13309655

# load GetOrganelle
get_organelle_from_reads.py -1 SRR13309655.1_1.fastq -2 SRR13309655.1_2.fastq -o SSR55 -R 15 -k 21,39,45,65,85,105,115,125 -F embplant_pt -t 20
get_organelle_from_reads.py -1 SRR13309655.1_1.fastq -2 SRR13309655.1_2.fastq -o SSR55_w075 -R 15 -k 21,39,45,65,85,105,115,125 -F embplant_pt -t 20 -w 0.75

# To run the GetOrganelle command, when the input data is too large, choose a fraction of the large dataset using seqtk module.
# example: J. validus SRR21976091_1.fastq and SRR21976091_2.fastq have about 87G, respectively.
module add seqtk
seqtk sample -s100 SRR21976091_1.fastq 0.05 >J_validus_s1.fastq
seqtk sample -s100 SRR21976091_2.fastq 0.05 >J_validus_s2.fastq
