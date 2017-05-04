# Thesis project: “Motif discovery in next-generation sequencing data”

Andrea Corneo, Politecnico di Milano

## Project context
This thesis project will focus on motif discovery, regarding in particular those
motifs bound by Transcription Factors (TF).

Discovering a “motif” means finding the sequence (8-12 base pairs (bp)
long) most likely to correspond to enriched or over-represented TF Binding
Site (TFBS) in thousands of DNA sequences.

A common approach to finding motifs is based on Chromatin Immunoprepicipation (ChIP) 
to extract DNA sequences (200 bp) bound by the TF. The
fragments are, then, sequenced to obtain small reads (50 bp). The latter are,
eventually, aligned on the whole genome to reconstruct the original sequence.
This technique is called ChIP-Seq and is followed by the Peak Calling phase, in
order to identify enriched DNA regions bound by the TF to be analyzed through
motif discovery tools.

## First part: motif discovery
The main novelty is based on the fact that motifs are already contained in the
original sequences, that is, before the Peak Calling phase, thus, an algorithm
should be able to find them directly in the not-preprocessed sequences. This
implies to work on millions or hundreds of millions sequences, instead of on
thousands of enriched regions coming from the Peak Calling phase.

We expect to find enriched motifs with high significance, since each instance
is sequenced more than once and thus appears more frequently then expected
inside the dataset.

The idea of performing motif discovery directly on the reads – if long enough,
otherwise on the mapped and extended ones – requires to use an efficient indexing method. While usually structures like suffix trees are sufficient, our
approach will be based on De Bruijn Graphs. This approach requires to adapt
commonly used motif discovery algorithms so that they will be able to work on
a much greater amount of shorter sequences.

## Second part: motif rules
This part aims at defining “recruitment rules” for the motifs found on DNA
regions bound by a TF, that is, how the TF interacts with other TFs. Regardless
of the method used for upstream motif discovery (the one proposed here or
others), we should expect the output to present several motifs identified as
enriched, each one with an associated score or significance measure.

Usually motif discovery ends here with a list of enriched motifs. Our approach will go further, aiming at finding correlations and anti-correlations of
enriched motifs in the regions bound by the TF.

This is, thus, a classification problem: we want to determine the rules and
partition DNA regions according to each rule. We will start to approach this
problem by using Hidden Markov Models (HMM).

## Datasets used
The analysis and validation of the methods developed will be performed on
data produced in the framework of the ENCODE Project1 . Using ChIP-Seq,
this project characterized hundreds of TFs in different cell lines; the quality
of these datasets constitutes an optimal benchmark for testing the proposed
methods.

In particular the validation of the methods developed will be performed on
datasets for TF NF-Y : for the first part raw data will be needed, while for the
definition of the rules it will start from Motif Discovery output – accordingly to
the fact that the HMM should be independent from the Motif Discovery step.

## Possible extension
An integration of this work with GenoMetric Query Language (GMQL) will be
evaluated later.

