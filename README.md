# metaOthello
## Installation

Download code from https://github.com/xa6xa6/metaOthello/

There are two code folders:

* Builder -- The code under "build" is used to generate the MetaOthello index.
* Classifier -- The code under "classifier" is used to perform taxonomic classification of sequencing reads using MetaOthello Index.

## Builder usage
  
### Ready-built indexes indexes 
Building indexing is very time-consuming (costs about 6 hours to build index for NCBI/refseq bacterial genome database).
Therefore, we provide ready-built indexes (NCBI/refseq bacterial genome database) for users to download:

* 20mer index: https://drive.google.com/open?id=0BxgO-FKbbXRIYWREa2NwejlVYUU
* 25mer index: https://drive.google.com/open?id=0BxgO-FKbbXRIY1pRaHJsYVg5dTQ
* 31mer index: https://drive.google.com/open?id=0BxgO-FKbbXRIa0Flc3Q4bWtycGM

### Building your own index
If you want to build an index with your own reference sequences, follow these steps. 

Jellyfish is used to prepare k-mer files for each reference sequence. Download *jellyfish* from: http://www.cbcb.umd.edu/software/jellyfish/
1. Produce a k-mer count file for your reference seqeuences. Command: 
```
    jellyfish count \
    –o <path_to_bacterial_rawKmerCountFile> \
    -m <Kmer_length> \
    -t <threads_num> \
    -s <bf_size> \
    -C <path_to_bacterial_referenceSeqFastaFile>
```
2. Dump k-mers to human-readable format. Command:
```
    jellyfish dump \
    –t –c \
    –o <path_to_bacterial_readableKmerCountFile> \
    <path_to_bacterial_rawKmerCountFile>
```
3. Generate taxonomy info file.
    Put all readable k-mer count files into the same directory 
    `<path_to_bacterial_reference_seq_Kmer_file_dir>`
    and rename them as `1.Kmer`, `2.Kmer`, `...`, `m.Kmer`, 
    and generate a taxonomy info file like: 
    https://drive.google.com/open?id=0BxgO-FKbbXRIZlV3ZzBBdlFpMTQ
    
    There are three columns for each taxonomic rank in the file: 
    the 1st column is a reissued id from `0` to `m-1`, 
    where `m` is the total taxon num in that taxonomic rank. 
    The 2nd column lists taxon ids, and the 3rd column lists taxon scientific names. 
    Each row represents a species and its associated taxonomy info.

4. Run `make build` under the directory `build`.

5. Build the index. Command:
```
    ./build \
    <bacterial_reference_seq_associated_taxonomy_info_file(generated in Step1.3)> \
    <path_to_bacterial_reference_seq_Kmer_file_dir> \
    <shared kmer file suffixes> \
    <Kmer_length> 6 \
    <path_to_bacterial_index> \
    <path_to_a_temp_dir_for_intermediate_files>
```

## Classifier usage

1. Run `make classifier` in the `classifier` directory. 
2. Perform taxonomic classification for each metagenomics sequencing reads. Command:
```
    ./classifier \
    <path_to_bacterial_index> \
    <path_to_output_results_dir> \
    <Kmer_length> \
    <threads_num> \
    <fa_or_fq> \
    <SE_or_PE> \
    <bacterial_speciesId2taxoInfo_file> \ 
    <NCBI_names_file> \
    <readFile_singleEnd or readFile_end1> \
    (<readFile_end2 if paired-end reads are provided>)
```

  `<bacterial_speciesId2taxoInfo_file>` can be downloaded from: 
  https://drive.google.com/open?id=0BxgO-FKbbXRIc3FkLVFvMlpVVGM    
    
  Each row represents a species and its associated taxon ids at each taxonomic rank:
  species, genus, family, order, class, and phylum. Assign `-1` if the taxon id is not available.

  `<NCBI_names_file>` can be downloaded from: 
  https://drive.google.com/open?id=0BxgO-FKbbXRIUFI2dHlBMXZhdTA

  **NOTE:** We will keep the following files updated with the latest NCBI/refseq bacterial genome databases:

 1. bacterial reference seq associated taxonomy info file,
 2. bacerial index (MetaOthello index for classification)
 3. bacterial speciesId2taxoInfo_file
 4. NCBI names file

  Also, we will release tools for generating all the above files (from NCBI/refseq bacterial genome databases) very soon.

# License

    Copyright (C) 2016-, University of Kentucky
    Please refer to LICENSE.TXT for the detailed 'License'.


