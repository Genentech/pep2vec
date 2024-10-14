# Pep2Vec: An Interpretable Model for Peptide-MHC Presentation Prediction and Contaminant Identification in Ligandome Datasets
## Abstract

As personalized cancer vaccines advance, precise modeling of antigen presentation by MHC class I (MHCI) and II (MHCII) is crucial. High quality training data is essential for models used in clinical settings. While existing deep learning models focus on prediction performance, they offer little interpretability to ensure confidence in the underlying data. Here, we introduce Pep2Vec, a modular, transformer-based model, trained on both MHC I and MHC II ligandome data, that transcends these limitations  by transforming each input sequence into an interpretable vector. This approach integrates source protein features and elucidates the source of the performance gains it enables, revealing regions that correlate with gene expression and protein-protein interactions. Furthermore, Pep2Vec’s peptide latent space manifests relationships between peptides of varying MHC class, MHC allotype, peptide lengths, and submotifs. This resolution enables the systematic identification of four major contaminant types which constitute 5.0% of our ligandome data. These innovations establish Pep2Vec as a performative model in MHC presentation prediction, offering not only heightened predictive accuracy but also reducing the recommendation of contaminant-like peptides. Pep2Vec addresses a critical need for the development of more precise and effective applications of peptide MHC models, such as for cancer vaccines and antibody deimmunization.


## Repository Information
This repo contains the binary for the "flagship" model from the pep2vec paper, as well as a dashboard for navigating the latent space of the public datasets used to train and test this model.  Note the embeddings of this latent space are generated from a model trained on the full dataset, whereas the binary here is trained on just the train portion of that dataset.   

## Visualization
See the 'viz' directory and the README.md there.

## Model Inference
### Installation
Simply copy the pep2vec.bin file to a local directory, it can be executed on any linux operating system.  
Note if checking the repository out make sure to install git lfs first to get the binary file.  
```bash
git lfs install
git clone ...
git lfs pull
```

### Usage
Inference is run with the following command:
```bash
./pep2vec.bin -–dataset /[full absolute path to repository]/test_input_mhc1.csv -–output_location [full absolute path to output folder]/test_output.parquet --mhctype mhc1 
```

Note that absolute paths are required for defining inputs and outputs.  The output is stored in a parquet format for efficiency, given there is a very large number of columns in the output file.

"mhctype" must be either mhc1 or mhc2 depending on the type of MHC being scored.

### Input Format
See input_example.csv for an example of the input format.  The input file should be a csv with the following columns:

    - n_flank
       - The 5 amino acid sequence immediately preceding the peptide, if unavailable use '*****'
    - peptide
       - The peptide sequence to be scored, should be between 9-25 amino acids long
    - c_flank
         - The 5 amino acid sequence immediately following the peptide, if unavailable use '*****'
    - allele
        - The MHC1 allele to be scored, should be in the format 'A*02:01__B*07:02' for a multi allele scoring those 2 alleles
        - See example file for more details
        - Leave empty if scoring MHC2
    - allotype
        - The MHC2 allotype to be scored, see the example file for the format for single allele or multi allele scoring
        - Leave empty if scoring MHC1
    - gene_ids
        - The ensemble gene id of the source protein
        - If multi maps to more than one source gene, use any gene_id with mapping
        - If unavailable leave empty, however this provided model has undefined performance on peptides without gene_ids

    

### Output Format
    - All input columns are retained in the output file
    - EL_pred
        - The predicted logit score of presentation likelihood
    - peptide_core
        - The predicted binding core of the peptide
    - latent_[0:1023]
        - The 1024 dimensional latent space representation of the peptide
    - attn_[0:50]
        - The 51 attention weights of the peptide
    - Additional column used for internal processing that can be ignored.
