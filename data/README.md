# Data
This folder contains raw data and materials which can be used to replicate findings.

## alignments

**extant_fluorescent_proteins_data.csv**
Table of over 70 documented natural fluorescent proteins, including name, species, sequence, original reference and, where possible, accession numbers and measured excitation/emission peaks.

**extant_fluorescent_proteins.fasta**
Fasta file of the above fluorescent proteins.

**extant_fluorescent_proteins_TCoffee_aligned.fasta**
Fasta file of the above proteins, structurally aligned (T-Coffee Expresso).

**avGFP_amacGFP_cgreGFP_ppluGFP2__protein_sequences.fasta**
Protein sequences of the four GFPs used in  the study.

**avGFP_amacGFP_cgreGFP_ppluGFP2__protein_sequences_TCoffee_aligned.fasta**
Structurally aligned sequences of the four GFPs used in the study (T-Coffee Expresso).


## cell_sorting

**`x`GFP__gate_border_values.txt**
Tables containing the absolute values between gate borders in the green channel during FACS.

### fcs_files
Folder containing recorded sorting data for amacGFP, cgreGFP, and ppluGFP2 genes, machines A and B (both FACS Aria III).


## final_datasets

**`x`GFP__barcodes_to_nucleotide_sequences.zip**
Fasta files linking primary (20N) barcodes to their corresponding full gene sequences, extracted from MiSeq data. Only barcodes passing quality control (see Methods: MiSeq data processing) are included.

**amacGFP_cgreGFP_ppluGFP2__final_nucleotide_genotypes_to_brightness.csv** and **amacGFP_cgreGFP_ppluGFP2__final_aminoacid_genotypes_to_brightness.csv**
Dataframes linking nucleotide or protein genotypes to their measured fluorescence level (see Methods). Mutations in genotypes are labeled in the format AiB, where A is the original wildtype state, B is the mutated state, and i is the position (counting starts from Methionine = 0). In the nucleotide dataset, 'n_replicates' refers to the combined number of distinct barcodes representing a genotype and machines it was measured on. In the amino acid dataset, 'n_replicates' refers to the number of synonymous nucleotide sequences measured for each protein sequence.

**avGFP__rf_nucleotide_genotypes_to_brightness.csv** and **avGFP__rf_aminoacid_genotypes_to_brightness.csv**
avGFP data from Sarkisyan et al 2016, reformatted to match the same indexing scheme as the new data (initial Methionine = 0).

**avGFP_amacGFP_cgreGFP_ppluGFP2__wt_and_synonymous_barcodes.csv**
Table of barcodes corresponding to sequences coding for wildtype proteins, for comparison of sequences with and without synonymous mutations.


## predictions

**experimentally_tested_predictions.csv**
Table containing sequences of ML-generated GFP genotypes, including their distance from the wildtype, predicted fluorescence according to two neural net models, and raw fluorescence values extracted from images of cells expressing the variants.


## protein_structure

**avGFP_amacGFP_cgreGFP_ppluGFP2__ddG_predictions.csv**
Predictions of delta-delta-G for oberserved single amino acid mutations in all proteins.

### PBD_structures 
Folder containing .pdb crystal structures for avGFP, amacGFP, cgreGFP, and ppluGFP2.

### residue_distance_matrices
Dataframes containing the minimum physical distance between pairs of residues inside the 3D GFP structures, in Angstroms. Row and column indices represent the residue position within the protein, starting from 0 for the initial methionine.

### secmals

**secmals_peak_weights.csv**
Calculated molecular weights of peaks detected using SEC-MALS on all wildtype proteins.

**secmals_raw.csv**
SEC-MALS data including retention volume, refractive index, ultraviolet, and right angle light scatter for all wildtype proteins.

### stability_measurements

**cd_melting_curves.csv**, **cd_spectra_200to260.csv**, and **circular_dichroism_raw_data.xlsx**
Circular dichroism raw data for melting curves (30-98C, 1C/minute ramp) or full-spectrum measurements (30C and 98C) of all wildtype proteins; monitored wavelengths are specified. In the .xlxs file, measurements for different proteins or different conditions are listed in separate tabs.

**dsc_data__av_amac_amacV14L_cgre_pplu.csv**
Differential scanning calorimetry data for all wildtype proteins.

**dsf_data__av_amac_amacV14L_cgre_pplu.csv**
Differential scanning fluorimetry data for all wildtype proteins. Unfolding/scattering curves and their derivatives are included.

**qpcr_melting_curves.csv**
Data from a melting curve (20-98C) conducted on all wildtype proteins.

**urea__absorbance_spectra.csv**
Table containing absorbance values of all genes in the range of 300 to 700nm, in 9M urea and PBS, at multiple consecutive time points (see Methods: Urea sensitivity assays). Blank control values are already subtracted.

**urea__fluorescence_spectra.csv**
Table containing fluorescence values of all genes in the range of 450 to 700nm (excitation at 420nm), in 9M urea and PBS, at multiple consecutive time points. Blank control values are already subtracted.


### unfiltered_barcodes_to_brightness

**`x`GFP__unfiltered_barcodes_to_brightness_machine_`X`.csv**
Dataframes containing the read distribution across gates of all primary-secondary barcode combinations, along with their fitted fitness values (see Methods). Read counts are normalized (see Methods), but data are not filtered according to cell count, number of replicates, etc.
