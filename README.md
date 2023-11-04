# mad_decoy

**Mad Decoy** is a Python script tool that allows the integration of protein q-values from multiple datasets. Mad decoy is used by quantms absolute quantification pipeline to integrate multiple datasets in the same view. This is needed because part of independent project FDR filtering when integrating multiple datasets the number of false positives could be inflated. 

In quantms absolute quant view (e.g https://quantms.org/ae/tissues?protein=LRBA_HUMAN), we apply the following QC rules for the integration: 

At the level of the dataset:
  
  1 - **1% FDR at PSM level**: with this filter, we guarantee that all the PSMs are high-quality.
  
  2 - **peptide AA > 7**: All peptide sequences identified must be higher than 7 amino acids.
  
  3 - **1% FDR at the protein level**: with this filter, we guarantee that all proteins reported are high-quality.
  
  4 - **2 unique peptides per protein**: In the previous step we statistically control the FDR at the protein level (1% FDR), extra filtering is applied and only proteins with 2 unique peptide sequences are reported.

**Note**: It is important to know that all datasets can be independently retrieved from quantms with only rules 1,2,3 applied for other types of analyses. 

However, when combining datasets, the number of false positives can increase really quickly. As many datasets get integrated more false positives accumulate. `mad_decoy` aims to provide multiple strategies and algorithms to statistically enable the integration of large-scale proteomics datasets while controlling the false positive accumulation. 

## Decoy simulation algorithm 

For the integration, of multiple datasets the first algorithm provided by **mad_decoy** tries to simulate for the list of proteins resulting from the integration datasets, a distribution of decoy proteins `fake decoys`. The `fake decoys` follow a q-value distribution similar to all the target proteins in the integration list. 

**Note**: This algorithm is needed because, in quantms tools like DIANN do not report back the decoy proteins quantified in the analysis, which does not allow the performance of more complex and accurate methods. 

```python
python add_fake_merge.py
```

The script must be executed in the folder where all the files containing the protein q-values are stored in `tsv` files. The structure for those files is the following: 

| protein_accessions  | protein_global_qvalue  |  is_decoy | condition | dataset_accession |
|---------------------|:----------------------:|----------:|----------:|------------------:|
| Q8NBF2              |  0.0003067249444061038 | 0         | platelets |   PXD000561       |
  
In the folder where the script will be run, all projects should be represented with the following structure. The script `add_fake_merge.py` will generate a list with all proteins in the integrated and the corresponding fake proteins and the adjusted q-value. 

## Contributors

Lukas Kall - Statistical Biotechnology at KTH

Yasset Perez-Riverol - EMBL-EBI 

