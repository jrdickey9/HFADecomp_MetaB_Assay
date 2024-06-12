#!/bin/bash

#file path to dir with folder of seqs
cd /Users/lab/Desktop/Jonathan_Dickey_LabMac/HFA/bacterial_prelim

###Qiime2 lines begin here
##Import
#qiime tools import \
#  --type 'SampleData[PairedEndSequencesWithQuality]' \
#  --input-path JRD_HFA_16S \
#  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
#  --output-path HFA-bac-paired-endR1.qza

##Visualize
#qiime demux summarize \
#  --i-data HFA-bac-paired-endR1.qza \
#  --o-visualization HFA-bac-paired-endR1.qzv

#using imported qza file to run the rest of the pipeline below. 
#FILTER, TRIM, AND DENOISE PAIRED DATA
#qiime dada2 denoise-paired \
#  --i-demultiplexed-seqs HFA-bac-paired-endR1.qza \
#  --p-trunc-len-r 132 \
#  --p-trunc-len-f 138 \
#  --p-max-ee-f 2 \
#  --p-max-ee-r 3 \
#  --p-pooling-method "independent" \
#  --p-chimera-method "consensus" \
#  --o-representative-sequences HFABac-paired-asv-sequencesR1.qza \
#  --o-table HFABac-paired-asv-tableR1.qza \
#  --o-denoising-stats HFABac-paired-denoising-statsR1.qza

##VISUALIZE PAIRED DATA
#qiime metadata tabulate \
#  --m-input-file HFABac-paired-denoising-statsR1.qza \
#  --o-visualization HFABac-paired-denoising-statsR1.qzv

#qiime feature-table tabulate-seqs \
#  --i-data HFABac-paired-asv-sequencesR1.qza \
#  --o-visualization HFABac-paired-asv-sequencesR1.qzv

#qiime feature-table summarize \
#  --i-table HFABac-paired-asv-tableR1.qza \
#  --m-sample-metadata-file updated_meta_bacterial_data.tsv \
#  --o-visualization HFABac-paired-asv-tableR1.qzv

qiime feature-table summarize \
  --i-table HFABac-paired-asv-tableR1.qza \
  --o-visualization HFABac-paired-asv-tableR1.qzv

#SILVA 138; TRAIN -- make sure its for 515-806
#qiime feature-classifier fit-classifier-naive-bayes \
#  --i-reference-reads silva-138-99-seqs-515-806.qza \
#  --i-reference-taxonomy silva-138-99-tax-515-806.qza \
#  --o-classifier silva-138-99-515F-806R-nb-classifier.qza

#CLASSIFY AGAINST SILVA w PAIRED DATA
#qiime feature-classifier classify-sklearn \
#  --i-classifier silva-138-99-515F-806R-nb-classifier.qza \
#  --i-reads HFABac-paired-asv-sequencesR1.qza \
#  --output-dir HFABac-paired-taxonomyFINAL

#VISUALIZE TAXONOMY of PAIRED DATA
#qiime metadata tabulate \
#  --m-input-file HFABac-paired-taxonomyFINAL/classification.qza \
#  --o-visualization HFABac-paired-taxonomyFINAL.qzv

#BUILD PHYLOGENY w PAIRED DATA 
#qiime phylogeny align-to-tree-mafft-fasttree \
#  --i-sequences HFABac-paired-asv-sequencesR1.qza \
#  --o-alignment HFABac-aligned-paired-seqsR1.qza \
#  --o-masked-alignment HFABac-masked-aligned-paired-seqsR1.qza \
#  --o-tree HFABac-unrooted-tree-paired-seqsR1.qza \
#  --o-rooted-tree HFABac-rooted-tree-paired-seqsR1.qza

#EXPORT ASV TABLE
#qiime tools export \
#  --input-path HFABac-paired-asv-tableR1.qza \
#  --output-path HFABac_paired_dadatableR1

#EXPORT TAXONOMY
#qiime tools export \
#  --input-path HFABac-paired-taxonomyFINAL/classification.qza \
#  --output-path HFABac_paired_taxaR9

#EXPORT PHYLOGENY
qiime tools export \
  --input-path HFABac-rooted-tree-paired-seqsR1.qza \
  --output-path HFABac-rooted-tree-paired-seqsR1







