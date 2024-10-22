# BCP-ALL_QE-TMT_Nat-Comm-2019

#### Re-analysis of data from childhood acute lymphoblastic leukemia study in Nat. Comm. April 2019

## Demonstrates analyzing isobaric labeling data without taking ratios

#### Phil Wilmarth, OHSU
#### _April 8, 2019_

---

### Paper Overview

The [paper](https://www.nature.com/articles/s41467-019-09469-3) is comparing samples from pediatric B-cell precursor acute lymphoblastic leukemia (BCP ALL) that were labeled with 10-plex TMT and analyzed on a Q-Exactive instrument. There were two leukemia conditions: high hyperdiploid (18 samples) and diploid/near-diploid ETV6/RUNX1-positive cases (9 samples). Samples were allocated in a balanced design with 6 and 3 samples per 10-plex, respectively. A pooled standard sample was also created and added to each plex (131N channel). Three 10-plex experiments were used to accommodate the 27 samples.

Digestion was via an eFASP protocol, peptide cleanup with an sp3 method, and TMT labeling was according to manufacturer's recommendations.

Each plex was first separated into 72 fractions by isoelectric point and each fraction run in a 60-min RP gradient. The QE was run in a top-10 mode at 70K resolution. MS2 HCD scans were acquired at 35K resolution.

> Yang, M., Vesterlund, M., Siavelis, I., Moura-Castro, L.H., Castor, A., Fioretos, T., Jafari, R., Lilljebjörn, H., Odom, D.T., Olsson, L. and Ravi, N., 2019. Proteogenomics and Hi-C reveal transcriptional dysregulation in high hyperdiploid childhood acute lymphoblastic leukemia. Nature communications, 10(1), p.1519.

---

### Re-analysis Overview

The [216 QE RAW files](https://www.ebi.ac.uk/pride/archive/projects/PXD010175) were downloaded and converted to MS2 file format using [MSConvert](http://proteowizard.sourceforge.net/) and the script from the [PAW pipeline](https://github.com/pwilmart/PAW_pipeline.git). There were about 1.4 million spectra per TMT 10-plex for a total of 4.3 million MS2 spectra. A basic TMT processing was done in the [PAW pipeline](https://github.com/pwilmart/PAW_pipeline.git) using:
- wider tolerance [Comet](http://comet-ms.sourceforge.net/) searches
- canonical UniProt human [reference proteome](https://github.com/pwilmart/fasta_utilities.git) (20764 sequences including contams)
- accurate mass conditioned score histograms
- target/decoy FDR method by peptide classes
- parsimonious protein inference with 2 peptides per protein
- homologous protein grouping
- quantification at the protein level using summed reporter ions
- [Jupyter notebook](https://jupyter.org/) statistical analysis in R using [edgeR](http://bioconductor.org/packages/release/bioc/html/edgeR.html)

---

### Notebook HTML files:

[balanced study averages](https://pwilmart.github.io/TMT_analysis_examples/Nat-Comm-2019_TMT_QE_averages.html) - slightly better IRS using plex averages

[single pooled standard](https://pwilmart.github.io/TMT_analysis_examples/Nat-Comm-2019_TMT_QE_pools.html) - single pooled internal standards have a little more uncertainty

---

### Proteomics analysis in publication

The starting data for the re-analysis is the data used in the paper and briefly described above. Each 10-plex was extensively fractionated and run on the fast scanning Q-Exactive platform. This is a deep proteome characterization no matter how you care to define that.

#### Key analysis points:
- sample processing, fractionation, and mass spec settings are in the publication
- Ensembl human protein database used
  - very complete so lots of peptide redundancy
- MS2 scans acquired at 35K resolution
  - maybe >50K is better to completely resolve N- and C-forms
- search engine was MSGF+
  - 10 ppm parent mass tolerance
  - other parameters appropriate for basic TMT experiment
- OpenMS IsobaricAnalyzer module for reporter ion extraction
- PSMs filtered to 1% FDR (no description of decoy database)
  - Percolator was used (maybe in an odd way?)
- gene symbol based protein grouping
- proteins filtered to a 1% FDR (minimum number of peptides per protein?)
- reporter ion ratios to pooled standard channel used for quantification
  - ratios formed at the PSM level (?)
  - only used PSMs unique to genes
  - median of ratios used for protein summary
  - only genes seen in all samples were quantified
- statistical DE analysis done with limma
  - log2 of median ratio
  - some sort of batch correction (batches can be factors in limma?)

#### Some open questions:

- was there any precursor interference filtering done? (there is none in PAW)
- were any reagent purity corrections used? (PAW does not have that option)
- what was the minimum peptide length cutoff? (PAW uses 7)
- how was ambiguity in protein to gene symbol mapping handled?
- were any sample loading type normalizations done before forming ratios?
- ratios to common channel should remove batch effects
  - were there any indications of batch effects? (IRS removes plex effects)

---

### PAW analysis key points

- wider tolerance Comet searches
  - parent ion tolerance was 1.25 Da
  - allows accurate mass to distinguish correct and incorrect matches
  - has sufficient number of scored sequences for deltaCN and other quantities
- canonical UniProt human reference proteome (20764 sequences including contams)
  - already curated to genes
  - isoforms are complicated to reliably identify
  - avoids some ambiguities in protein inference
- accurate mass conditioned score histograms
  - improves sensitivity
  - catches deamidated Asn and first isotopic triggers
  - catches MS2 scans without MS1 information sufficient for an accurate mass call  
- target/decoy FDR method by separate peptide classes
  - accurate mass window
  - charge state (2+, 3+, and 4+ only)
  - modification status (unmodified and oxidized Met)
- parsimonious protein inference with 2 peptides per protein
  - identical peptide sets grouped together
  - formal peptide subsets removed
  - 2 peptides per protein per sample required
- homologous protein grouping
  - extends parsimony to include *nearly identical* and *nearly a subset*
  - combines proteins with mostly common shared peptides into families
- quantification at the protein level using summed reporter ions
  - only PSMs from peptides unique to final protein list are used
  - reporter ions summed into protein totals per plex before any normalizations
- [Internal reference scaling](https://github.com/pwilmart/IRS_validation.git) method used to match intensities between plexes
  - internal reference channels can be used (better if more than one)
  - plex-wide average can be [used](https://github.com/pwilmart/IRS_normalization.git) in balanced study designs
- statistical analysis in R using edgeR
  - supports common study designs
  - has nice built-in normalization (TMM)
  - uses variance sharing between genes/proteins for low-replicate number studies
  - visualizations (MA plots and MDS clustering)
  - multiple testing corrections

---

### Some interesting things

PSM identifications (1% FDR): PAW pipeline seemed to have 80% more identified PSMs than the numbers of identified PSMs listed in the paper's Supplemental file. This is a large factor and it seems too good to be true.

|Plex|Paper Supp. File|PSM file from PRIDE|PAW Analysis|
|----|-----------|------------|------|
|A|237K|461K|409K|
|B|217K|447K|404K|
|C|191K|386K|352K|
|Total|645K|1293K|1166K|

The numbers of identified PSMs in the Supplemental file of the paper (645K) did not look right. There was a PSM results table in the PRIDE archive along with the RAW files. That was processed to count PSMs and those number were added to the table above. The 1.29 million PSM number looks better. MSGF+ is considered a very good performing search engine and Percolator is a good classifier. Now, the correct result is that the PAW PSM identification number was less than the paper's number (90% of their number).

It is always tricky to know if counting of things for comparisons are being done under the exact same conditions. There are lots of little filters (charge states, peptide lengths, enzymatic status, etc.) that can be done at many points in pipelines. The number of unique peptides reported in the paper was 175K. The PAW number at the step where reporter ions are summed to protein totals was 172K.

Protein numbers are also an apples-to-oranges situation. There were different protein databases and summarization methods used. These factors affect the union of the protein IDs more strongly that they affect the intersection of identifications across the three plexes. Observation in all three plexes is a noise filter that removes some of the lowest abundance proteins. That will tighten up the protein FDRs and improve comparisons. The paper had 8480 quantifiable proteins seen in all samples. There were 8756 quantifiable proteins with the PAW processing.

The numbers of differential candidates were also a little different. The paper has 1286 up regulated and 1127 down regulated proteins at 5% FDR. The PAW edgeR results using reporter ions in their natural intensity scale rather than as transformed ratios, had 962 and 970 at a 10% FDR. Statistical testing p-values, multiple testing correction methods, and the numbers of DE candidates at arbitrary cutoffs are much more dynamic that people realize. To me, the agreement looks pretty good. These are one-protein-at-a-time DE assessments. Most biological processes involve more than one protein. Biological interpretation of differential expression needs to be done with sets of proteins. We have very few, very primitive tools for this currently.

The single pooled standard was noisier that the plex-wide averages for the IRS procedure. Protein coefficients of variation (CVs) were better and the statistical testing produced a small increase in candidates when the averages were used instead of the single common channel.

---

## Summary

- This is a deep dive into the cancer B-cell proteome with almost 50% of the reference proteome detected
- The TMT data from the QE for this experiment looked okay
  - the SPS MS3 method was developed to reduce interference and compression in MS2 reporter ions
  - data compressed might be compensated by sample groups that have large expression differences
    - message expression range was greater than protein expression range in paper
- visualizations in R give a much more convincing view of the data
- notebooks that encapsulate the full linear analysis of data can tell a full story
  - traditional Supplemental files do not communicate sufficient details and information
- the experiment and data are quite nice
