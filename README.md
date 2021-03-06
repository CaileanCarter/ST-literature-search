[![Total alerts](https://img.shields.io/lgtm/alerts/g/CaileanCarter/ST-literature-search.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/CaileanCarter/ST-literature-search/alerts/)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/CaileanCarter/ST-literature-search.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/CaileanCarter/ST-literature-search/context:python)

# ST-literature-search

## What is it?
ST Literature Search is a command-line tool I use to quickly find literature relating to any given <i>E. coli</i> sequence type. 

## Requirements
- Python 3.6 or later
- Pandas
- Numpy
- BioPython

## How does it work?
If you search "Escherichia coli sequence type" into [PubMed](https://pubmed.ncbi.nlm.nih.gov/) you will find >25,000 results. Use the save options to save all the results as a CSV file (PubMed will only allow you to save the first 10,000 results). 
<br>
Using this tool, run the following command: 
```
python litsearch.py new [file_path] [email address]
```
<b>Tip:</b> you will need to provide your email address. As this tool uses the BioPython Entrez API, an email address is required to ensure PubMed does not see you as a bot. If no email is provided, you will be prompted to provide one. This tool does not store your email address or give it to any third-party other than the BioPython Entrez API.

<br>

This filters out all literature results which do not mention a sequence type in the title (other selections methods may be considered in the future). 

### Find papers
You can identify papers mentioning a given sequence type using the `find` argument:
```
python litsearch.py find [file_path] [ST###]
```
If your file path is simply "escherichi-set.csv" which was the default name given by PubMed during download, then you won't need to add the file_path argument
<br>
Example:
```
>>> python litsearch.py find ST131
PMID                                                  Title                                            Authors                                           Citation                           DOI            ST  Times cited
18077311  Intercontinental emergence of Escherichia coli...  Nicolas-Chanoine MH, Blanco J, Leflon-Guibout ...  J Antimicrob Chemother. 2008 Feb;61(2):273-81....            10.1093/jac/dkm464         ST131          288
21081548  Escherichia coli O25b-ST131: a pandemic, multi...               Rogers BA, Sidjabat HE, Paterson DL.  J Antimicrob Chemother. 2011 Jan;66(1):1-14. d...            10.1093/jac/dkq415         ST131          245
20572763  Escherichia coli sequence type ST131 as the ma...  Johnson JR, Johnston B, Clabots C, Kuskowski M...  Clin Infect Dis. 2010 Aug 1;51(3):286-94. doi:...                10.1086/653932         ST131          213
24345742  The epidemic of extended-spectrum-β-lactamase-...  Price LB, Johnson JR, Aziz M, Clabots C, Johns...  mBio. 2013 Dec 17;4(6):e00377-13. doi: 10.1128...         10.1128/mBio.00377-13         ST131          171
19474064  Rapid detection of the O25b-ST131 clone of Esc...  Clermont O, Dhanji H, Upton M, Gibreel T, Fox ...  J Antimicrob Chemother. 2009 Aug;64(2):274-7. ...            10.1093/jac/dkp194         ST131          109
27006459  Evolutionary History of the Global Emergence o...  Stoesser N, Sheppard AE, Pankhurst L, De Maio ...  mBio. 2016 Mar 22;7(2):e02162. doi: 10.1128/mB...         10.1128/mBio.02162-15         ST131          109
19687243  Complete nucleotide sequences of plasmids pEK2...  Woodford N, Carattoli A, Karisik E, Underwood ...  Antimicrob Agents Chemother. 2009 Oct;53(10):4...          10.1128/AAC.00688-09         ST131          104
22053197  Insights into a multidrug resistant Escherichi...  Totsika M, Beatson SA, Sarkar S, Phan MD, Pett...  PLoS One. 2011;6(10):e26578. doi: 10.1371/jour...  10.1371/journal.pone.0026578         ST131           92
23926176  Escherichia coli sequence type 131 (ST131) sub...  Colpan A, Johnston B, Porter S, Clabots C, Anw...  Clin Infect Dis. 2013 Nov;57(9):1256-65. doi: ...            10.1093/cid/cit503  ST131, ST131           87

[10 rows x 12 columns]
```

### Get paper summary
You can get the summary page for any paper on PubMed using the argument `summary` and providing the PMID (see above):
```
python litsearch.py summary [file_path] [PMID]
```
Example:
```
python litsearch.py summary 23926176

1. Clin Infect Dis. 2013 Nov;57(9):1256-65. doi: 10.1093/cid/cit503. Epub 2013 Aug
6.

Escherichia coli sequence type 131 (ST131) subclone H30 as an emergent
multidrug-resistant pathogen among US veterans.

Colpan A(1), Johnston B, Porter S, Clabots C, Anway R, Thao L, Kuskowski MA,
Tchesnokova V, Sokurenko EV, Johnson JR; VICTORY (Veterans Influence of Clonal
Types on Resistance: Year 2011) Investigators.

Collaborators: Allen BL, Baracco GJ, Bedimo R, Bessesen M, Bonomo RA, Brecher SM,
Brown ST, Castellino L, Desai AS, Fernau F, Fisher MA, Fleckenstein J, Fleming
CS, Fries NJ, Kan VL, Kauffman CA, Klutts S, Ohl M, Russo T, Swiatlo A, Swiatlo
E.

Author information:
(1)University of Minnesota, Minneapolis.

Comment in
    Clin Infect Dis. 2013 Nov;57(9):1266-9.

BACKGROUND: Escherichia coli sequence type 131 (ST131), typically
fluoroquinolone-resistant (FQ-R) and/or extended-spectrum β-lactamase
(ESBL)-producing, has emerged globally. We assessed its prevalence and
characteristics among US veterans.
METHODS: In 2011, 595 de-identified E. coli clinical isolates were collected
systematically within 3 resistance groups (FQ-susceptible [FQ-S], FQ-R, and
ESBL-producing) from 24 nationally distributed Veterans Affairs Medical Centers
(VAMCs). ST131 and its H30 subclone were detected by polymerase chain reaction
and compared with other E. coli for molecular traits, source, and resistance
profiles.
RESULTS: ST131 accounted for 78% (184/236) of FQ-R and 64.2% (79/123) of
ESBL-producing isolates, but only 7.2% (17/236) of FQ-S isolates (P < .001). The
H30 subclone accounted for ≥95% of FQ-R and ESBL-producing, but only 12.5% of
FQ-S, ST131 isolates (P < .001). By back-calculation, 28% of VAMC E. coli
isolates nationally represented ST131. Overall, ST131 varied minimally in
prevalence by specimen type, inpatient/outpatient source, or locale; was the most
prevalent ST, followed distantly by ST95 and ST12 (13% each); and accounted for
≥40% (β-lactams), >50% (trimethoprim-sulfamethoxazole , multidrug), or >70%
(ciprofloxacin, gentamicin) of total antimicrobial resistance. FQ-R and
ESBL-producing ST131 isolates had higher virulence scores than corresponding
non-ST131 isolates. ST131 pulsotypes overlapped extensively among VAMCs.
CONCLUSIONS: Among US veterans, ST131, primarily its H30 subclone, accounts for
most antimicrobial-resistant E. coli and is the dominant E. coli strain overall.
Possible contributors include multidrug resistance, extensive virulence gene
content, and ongoing transmission. Focused attention to ST131, especially its H30
subclone, could reduce infection-related morbidity, mortality, and costs among
veterans.

DOI: 10.1093/cid/cit503
PMCID: PMC3792724
PMID: 23926176  [Indexed for MEDLINE]
```

Link to paper via DOI will be functional.