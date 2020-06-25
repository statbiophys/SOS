#### Simple Olga Sonia (SOS)

SOS evaluates the probability of generation (P<sub>gen</sub>) and probability in the periphery (P<sub>post</sub>) of specific T and B cell receptor sequences of humans and mice.

SOS also generates synthetic generated and peripheral repertoires of naive cells.

SOS is built on top of the [IGOR](https://github.com/qmarcou/IGoR)\[2\], [OLGA](https://github.com/zsethna/OLGA)\[3\] and [SONIA](https://github.com/statbiophys/SONIA/)\[4,5\] software.

P<sub>gen</sub> is conditioned on productivity of the sequence. It only takes into account VDJ recombination in the absence of selection. P<sub>post</sub> takes selection into account, and should better reflect frequency of occurrence in the periphery. However note that it is insensitive to HLA type. The models were learned on a possibly biased panel of individuals and may depend on the sequencing technologies used. For more accurate dataset-specific comparisons it is advised to infer new models using IGOR and SONIA.

<br />


Evaluate command:
`sonia-evaluate --humanTRB CASSTGNYGAFF --v_mask TRBV9 --j_mask TRBJ-1`

<br />

Generate command:
`sonia-generate --humanTRB -n 1000 --pgen (or --ppost)`

<br />

To cite this tool:
\[1\] Isacchini et al (2020) SOS: A web server for SONIA and OLGA, in preparation.

##### References:
- [\[2\] Marcou Q, Mora T, Walczak AM (2018) High-throughput immune repertoire analysis with IGoR. Nat Commun 9, 561 (2018)](https://www.nature.com/articles/s41467-018-02832-w)
- [\[3\] Sethna Z, Elhanati Y, Callan CG, Walczak AM, Mora T(2019) OLGA: fast computation of generation probabilities of B- and T-cell receptor amino acid sequences and motifs. Bioinformatics, Volume 35, Issue 17, 1 September 2019, Pages 2974–2981](https://academic.oup.com/bioinformatics/article/35/17/2974/5292315)
- [\[4\] Sethna Z, Isacchini G, Dupic T, Mora T, Walczak AM, Elhanati Y, Population variability in the generation and thymic selection of T-cell repertoires, (2020) bioRxiv, https://doi.org/10.1101/2020.01.08.899682](https://www.biorxiv.org/content/10.1101/2020.01.08.899682v1)
- [\[5\] Isacchini G,Sethna Z, Elhanati Y ,Nourmohammad A, Mora T, Walczak AM, On generative models of T-cell receptor sequences,(2019) bioRxiv, https://doi.org/10.1101/857722](https://www.biorxiv.org/content/10.1101/857722v2)
