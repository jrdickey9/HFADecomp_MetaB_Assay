#!/bin/bash

cd /Users/lab/Desktop/Jonathan_Dickey_LabMac/HFA/fungal_prelim 

qiime dada2 denoise-single \
  --i-demultiplexed-seqs HFAFung-sequencesR13.qza \
  --p-trunc-len 235 \
  --p-max-ee 2 \
  --o-representative-sequences HFAFung-forward-asv-sequencesR1.qza \
  --o-table HFAFung-forward-asv-tableR1.qza \
  --o-denoising-stats HFAFung-forward-denoising-statsR1.qza

qiime metadata tabulate \
  --m-input-file HFAFung-forward-denoising-statsR1.qza \
  --o-visualization HFAFung-forward-denoising-statsR1.qzv

qiime feature-table tabulate-seqs \
  --i-data HFAFung-forward-asv-sequencesR1.qza \
  --o-visualization HFAFung-forward-asv-sequencesR1.qzv

qiime feature-table summarize \
  --i-table HFAFung-forward-asv-tableR1.qza \
  --m-sample-metadata-file updated_meta_fungal_data.tsv \
  --o-visualization HFAFung-forward-asv-tableR1.qzv

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path sh_refs_qiime_ver9_dynamic_25.07.2023.fasta \
  --output-path unite.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path sh_taxonomy_qiime_ver9_dynamic_25.07.2023.txt \
  --output-path unite-taxonomy.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads unite.qza \
  --i-reference-taxonomy unite-taxonomy.qza \
  --o-classifier HFAFung-classifier.qza

qiime feature-classifier classify-sklearn \
  --i-classifier HFAFung-classifier.qza \
  --i-reads HFAFung-forward-asv-sequencesR1.qza \
  --o-classification HFAFung-forward-taxonomyR1.qza

qiime metadata tabulate \
  --m-input-file HFAFung-forward-taxonomyR1.qza \
  --o-visualization HFAFung-forward-taxonomyR1.qzv

qiime tools export \
  --input-path HFAFung-forward-asv-tableR1.qza \
  --output-path HFAFung_forward_dadatableR1

qiime tools export \
  --input-path HFAFung-forward-taxonomyR1.qza \
  --output-path HFAFung_forward_taxaR1


