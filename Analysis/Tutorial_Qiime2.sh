#!/bin/bash

cd /Users/lab/Documents/Cyrus_Daruwala/Microplastics2/

#IMPORT DATA;
#qiime tools import \
#  --type 'SampleData[PairedEndSequencesWithQuality]' \
#  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
#  --input-path seqs/ \
#  --output-path Final_Qiime2_outputs/mp2-sequences.qza

cd /Users/lab/Documents/Cyrus_Daruwala/Microplastics2/Final_Qiime2_outputs

#VISUALIZE DATA IMPORT;
#qiime demux summarize \
# --i-data mp2-sequences.qza \
# --o-visualization mp2-sequences.qzv

#PAIRED END DATA (Forward and Reverse) and DADA2;
#qiime dada2 denoise-paired \
#  --i-demultiplexed-seqs mp2-sequences.qza \
#  --p-trunc-len-r 140 \
#  --p-trunc-len-f 145 \
#  --p-max-ee-f 2 \
#  --p-max-ee-r 2 \
#  --p-pooling-method "pseudo" \
#  --p-chimera-method "consensus" \
#  --o-representative-sequences mp2-paired-asv-sequencesFINAL.qza \
#  --o-table mp2-paired-asv-tableFINAL.qza \
#  --o-denoising-stats mp2-paired-denoising-statsFINAL.qza

#VISUALIZE PAIRED END DATA;
#qiime metadata tabulate \
#  --m-input-file mp2-paired-denoising-statsFINAL.qza \
#  --o-visualization mp2-paired-denoising-statsFINAL.qzv

#qiime feature-table tabulate-seqs \
#  --i-data mp2-paired-asv-sequencesFINAL.qza \
#  --o-visualization mp2-paired-asv-sequencesFINAL.qzv

#qiime feature-table summarize \
#  --i-table mp2-paired-asv-tableFINAL.qza \
#  --m-sample-metadata-file 16S_mp2_metadata.tsv \
#  --o-visualization mp2-paired-asv-tableFINAL.qzv

#SILVA 138; TRAIN -- make sure its for 515-806 #NEVER RUN THIS CHUNK START ON CLASSIFY AGAINST SILVA, #PAIRED sdfghjk
#qiime feature-classifier fit-classifier-naive-bayes \
#  --i-reference-reads silva-138-99-seqs-515-806.qza \
#  --i-reference-taxonomy silva-138-99-tax-515-806.qza \
#  --o-classifier silva-138-99-515F-806R-nb-classifierFINAL.qza

#CLASSIFY AGAINST SILVA, PAIRED
qiime feature-classifier classify-sklearn \
 --i-classifier silva-138-99-515F-806R-nb-classifierFINAL.qza \
 --i-reads mp2-paired-asv-sequencesFINAL.qza \
 --output-dir mp2-paired-taxonomyFINAL2

#VISUALIZE TAXONOMY PAIRED
qiime metadata tabulate \
  --m-input-file mp2-paired-taxonomyFINAL2/classification.qza \
  --o-visualization mp2-paired-taxonomyFINAL2.qzv

#BUILD PHYLOGENY PAIRED 
#qiime phylogeny align-to-tree-mafft-fasttree \
#  --i-sequences mp2-paired-asv-sequencesFINAL.qza \
#  --o-alignment mp2-aligned-paired-seqsFINAL.qza \
#  --o-masked-alignment mp2-masked-aligned-paired-seqsFINAL.qza \
#  --o-tree mp2-unrooted-tree-paired-seqsFINAL.qza \
#  --o-rooted-tree mp2-rooted-tree-paired-seqsFINAL.qza

#EXPORT PAIRED
#qiime tools export \
#  --input-path mp2-paired-asv-tableFINAL.qza \
#  --output-path mp2_paired_dadatableFINAL

qiime tools export \
  --input-path mp2-paired-taxonomyFINAL2/classification.qza \
  --output-path mp2_paired_taxaFINAL2

#EXPORT PHYLOGENY
#qiime tools export \
#  --input-path mp2-rooted-tree-paired-seqsFINAL.qza \
#  --output-path mp2-rooted-tree-paired-seqsFINAL

