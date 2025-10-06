# 🧬 Whole Genome Variant Analysis — *Pseudomonas aeruginosa* (SRR34663677)

![Bioinformatics](https://img.shields.io/badge/Field-Bioinformatics-blue)
![NGS](https://img.shields.io/badge/Analysis-WGS%20Variant%20Calling-orange)
![Tools](https://img.shields.io/badge/Tools-FastQC%20|%20BWA%20|%20GATK%20|%20SnpEff-green)
![Status](https://img.shields.io/badge/Status-Completed-success)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

> **Comprehensive Whole Genome Sequencing (WGS) variant analysis** of *Pseudomonas aeruginosa* clinical isolate (**SRR34663677**) using a reference-based workflow — from **raw read QC** to **biological interpretation and antimicrobial resistance (AMR) insights**.

---

## 🧭 Project Overview

This repository presents a **complete bacterial variant analysis pipeline**, covering:

- Retrieval of sequencing data from **NCBI SRA**
- **Quality Control** and trimming of raw FASTQ files
- **Read alignment** to a reference genome (PAO1)
- **Variant calling** using GATK HaplotypeCaller
- **Annotation** with SnpEff
- **Biological interpretation** and comparative analysis using NCBI Pathogen Detection Browser

---

**Pipeline:**  
`Raw Reads → Quality Control → Alignment → Variant Calling → Annotation → Interpretation`

---

## 📂 Dataset Information

| Attribute | Description |
|------------|--------------|
| **Organism** | *Pseudomonas aeruginosa* (clinical isolate) |
| **Accession ID** | SRR34663677 |
| **Sequencing Type** | Illumina paired-end (151 bp) |
| **Reference Genome** | *P. aeruginosa* PAO1 (RefSeq: GCF_000006765.1) |
| **Genome Size** | 6.26 Mb |
| **GC Content** | 66% |

---

## ⚙️ Workflow Steps
 🔹 1. Data Retrieval

Raw reads downloaded using NCBI SRA Toolkit:

```bash
fasterq-dump SRR34663677
Converted .sra to paired FASTQ files:

SRR34663677_1.fastq.gz
SRR34663677_2.fastq.gz

🔹 2. Quality Control & Trimming

Tools: FastQC, Trim Galore

Pre-trimming stats:

Total reads: 1,436,195

GC: 66%

Sequence length: 151 bp

No poor-quality reads detected

Post-trimming stats:

Reads retained: 1,413,659 (~98.4%)

Quality improved (Phred > 30)

Adapters removed

GC bias consistent with P. aeruginosa biology

Interpretation:
Trimmed data is clean and suitable for mapping and variant detection.

🔹 3. Genome Alignment

Tool: BWA-MEM, SAMtools

bwa index GCF_000006765.1_ASM676v1_genomic.fna
bwa mem -t 8 GCF_000006765.1_ASM676v1_genomic.fna \
SRR34663677_1_val_1.fq.gz SRR34663677_2_val_2.fq.gz > SRR34663677_aligned.sam


SAM → BAM → Sorted → Mark Duplicates:

samtools view -Sb SRR34663677_aligned.sam | samtools sort -o SRR34663677_sorted.bam
samtools markdup SRR34663677_sorted.bam SRR34663677_markdup.bam
samtools index SRR34663677_markdup.bam


Results Summary:

Metric	Value
Total Reads	2,778,805
Mapped Reads	86.4%
Properly Paired	85.9%
Avg. Coverage	~56×
Avg. MAPQ	35.1
GC (mapped)	66%

Interpretation:
Excellent mapping quality. Unmapped reads likely correspond to strain-specific genes or mobile genetic elements.

🔹 4. Variant Calling (GATK)

Tool: GATK v4.2 HaplotypeCaller

gatk HaplotypeCaller \
  -R GCF_000006765.1_ASM676v1_genomic.fna \
  -I SRR34663677_markdup.bam \
  -O SRR34663677_raw_variants.g.vcf.gz \
  --sample-ploidy 1 \
  -ERC GVCF


Filtering:

gatk VariantFiltration \
  -V SRR34663677_raw_variants.g.vcf.gz \
  -O SRR34663677_filtered.vcf.gz \
  --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
  --filter-name "BacterialHardFilter"


Variant Statistics:

High Ts/Tv ratio (3.15)

Mean depth: 60–70×

Majority variants = PASS

Allele Frequency (AF) = 1.00 (haploid isolate)

🔹 5. Variant Annotation (SnpEff)

Tool: SnpEff

java -jar snpEff.jar build -gff3 Pseudomonas_aeruginosa_PA01
java -jar snpEff.jar ann Pseudomonas_aeruginosa_PA01 SRR34663677_filtered.vcf.gz > SRR34663677_annotated.vcf


Annotation Summary:

Impact	Type	% Frequency
MODIFIER	Intergenic / upstream variants	80%
LOW	Synonymous	18%
MODERATE	Missense	1.7%
HIGH	Stop-gained / frameshift	0.1%

Example Genes Affected:

dnaA, recF, gyrB – DNA replication & repair

mutS, mexR, ampC – Antibiotic resistance mechanisms

Codon/Amino Acid Substitutions:

Common silent: CTG→CTG, GCG→GCG

Missense: ATC→GTC (I→V), CGC→CTC (R→L)

Nonsense: CAG→TAG (stop gained)

🔹 6. Functional & Biological Interpretation

Findings:

Synonymous variants (~78%): Neutral but can affect codon bias.

Missense variants (~1.7%): May alter protein function or stability.

Frameshift/Nonsense (0.1%): Potential truncations affecting regulatory proteins.

Regulatory variants: Could influence promoter strength or gene expression.

Functional impact:
Variants in ampC, mexR, and gyrB suggest possible antimicrobial resistance (AMR) modulation.

🔹 7. Comparative Analysis — NCBI Pathogen Detection

The annotated genome was compared with the NCBI Pathogen Detection Browser
.

Cluster ID	Total Isolates	Source	SNP Difference	AMR Genes
PDS000095640.50	597	Clinical & Environmental	0	aadA11, aph(3’)-Ib, gyrA_T83I, parC_S87W
PDS000075445.200	268	Clinical	<5	–
PDS000076787.33	116	Environmental	<10	–

Interpretation:

The isolate clusters with widely distributed global strains.

Presence of fluoroquinolone resistance mutations (gyrA, parC).

Indicates potential for clinical relevance and surveillance tracking.

🧠 Key Insights

Complete WGS-based variant pipeline validated on bacterial isolate.

Reliable identification of SNPs, indels, and functionally significant mutations.

Integration with public repositories enhances reproducibility and epidemiological context.

Supports AMR prediction, evolutionary studies, and clinical genomics research.

🧰 Tools & Technologies
Category	Tool
QC & Trimming	FastQC, Trim Galore
Alignment	BWA-MEM, SAMtools
Variant Calling	GATK (HaplotypeCaller)
Annotation	SnpEff
Comparative Analysis	NCBI Pathogen Detection Browser
Environment	Linux, Bash CLI
🧾 Keywords

Bioinformatics Whole Genome Sequencing Variant Calling Pseudomonas aeruginosa
GATK FastQC BWA-MEM SnpEff Antimicrobial Resistance
Annotation Comparative Genomics Bacterial Genomics

👩‍💻 Author

Debopriya
Bioinformatics & Genomics Enthusiast

🔗 https://github.com/DEBOPRIYA2320

📧 debopriya0920@gmail.com

📜 License

This project is licensed under the MIT License — free for research and educational use.

⭐ Acknowledgements

NCBI SRA

GATK Toolkit

SnpEff

FastQC

BWA

NCBI Pathogen Detection Browser

Special thanks to the open-source bioinformatics community for their tools and data accessibility.

⭐ If you find this project useful, please give it a star and share it with fellow researchers!
