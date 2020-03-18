# Simple Olga & Sonia (SOS)

Source code of the SOS web server at: https://sites.google.com/view/statbiophysens/sos 

Command to run locally on the terminal: python sos.py 

== Description ==

SOS evaluates the probability of generation (Pgen) and probability in the periphery (Ppost) of specific T and B cell receptor sequences of humans and mice.

SOS also generates synthetic generated and peripheral repertoires of naive cells.

SOS is built on top of the IGOR [1], OLGA [2] and SONIA [3,4] software.

Pgen is conditioned on productivity of the sequence. It only takes into account VDJ recombination in the absence of selection. Ppost takes selection into account, and should better reflect frequency of occurrence in the periphery. However note that it is insensitive to HLA type. The models were learned on a possibly biased panel of individuals and may depend on the sequencing technologies used. For more accurate dataset-specific comparisons it is advised to infer new models using IGOR and SONIA.

== To cite this tool ==

Isacchini et al (2020) A web server for SONIA and OLGA, in preparation.  

== References ==

[1] Marcou Q, Mora T, Walczak AM (2018) High-throughput immune repertoire analysis with IGoR.Nature Communication s9:561

[2] Sethna Z, Elhanati Y, Callan CG, Walczak AM, Mora T(2019) OLGA: fast computation of generation probabilities of B- and T-cell receptor amino acid sequences and motifs.Bioinformaticsbtz035

[3] Sethna Z, Isacchini G, Dupic T, Mora T, Walczak AM, Elhanati Y, Population variability in the generation and thymic selection of T-cell repertoires

[4] Isacchini G,Sethna Z, Elhanati Y ,Nourmohammad A, Mora T, Walczak AM, On generative models of T-cell receptor sequences


== Python dependencies  ==

- olga

- sonia

- dash

- das_daq

- pandas 

- numpy