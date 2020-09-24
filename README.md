# Simple Olga & Sonia (SOS)
## Description

SOS evaluates the probability of generation (Pgen) and probability in the periphery (Ppost) of specific T and B cell receptor sequences of humans and mice. SOS also generates synthetic generated and peripheral repertoires of naive cells. SOS is built on top of the IGOR [1], OLGA [2] and SONIA [3,4] software.

Pgen is conditioned on productivity of the sequence. It only takes into account VDJ recombination in the absence of selection. Ppost takes selection into account, and should better reflect frequency of occurrence in the periphery. However note that it is insensitive to HLA type. The models were learned on a possibly biased panel of individuals and may depend on the sequencing technologies used. For more accurate dataset-specific comparisons it is advised to infer new models using IGOR and SONIA.

### Run Online:
Source code of the SOS web server at: https://sites.google.com/view/statbiophysens/sos 


### Run Locally:
Install python dependencies listed below.
Change the name of the source of the local directory listed in the `utils.py` file. The path source should be named as the local directory contianing the SOS-Master folder on your computer. Ex:
`local_directory='/Users/Johny_Appleseed/SOS-master/'`

Command to run locally on the terminal: `python sos.py`

After running this command the system should return a local address which can be entered into the search bar of your internet browser to run the dash app.
Ex: `Dash is running on http://127.0.0.1:8050/`

## Python dependencies

- olga

- sonia

- dash

- dash_daq

- pandas 

- numpy==1.18.0

## To cite this tool

Isacchini et al (2020) SOS: online probability estimation and generation of T-and B-cell receptors,  Bioinformatics, btaa574, https://doi.org/10.1093/bioinformatics/btaa574

## References

[1] Marcou Q, Mora T, Walczak AM (2018) High-throughput immune repertoire analysis with IGoR.Nature Communication s9:561

[2] Sethna Z, Elhanati Y, Callan CG, Walczak AM, Mora T(2019) OLGA: fast computation of generation probabilities of B- and T-cell receptor amino acid sequences and motifs.Bioinformaticsbtz035

[3] Sethna Z, Isacchini G, Dupic T, Mora T, Walczak AM, Elhanati Y, Population variability in the generation and thymic selection of T-cell repertoires

[4] Isacchini G,Sethna Z, Elhanati Y ,Nourmohammad A, Mora T, Walczak AM, On generative models of T-cell receptor sequences
