# RoLiM - Robust Linear Motif Deconvolution
Detection of statistically enriched patterns in protein sequence data sets.

## 🔔 Important Notice Regarding System Updates and Limitations

Due to recent changes in our server infrastructure, we want to make you aware of some current limitations in our email notification system:

- System memory constraints may occasionally prevent email notifications from being sent, even when your analysis is completed successfully
- Results files exceeding 15MB cannot be sent via email, though they remain accessible through our platform

We have identified additional issues with timeout and error processes when using non-human or isoform FASTA files. If your analysis requires custom FASTA files, this functionality is temporarily unavailable. 

We understand these limitations may impact your workflow, and we're actively working on solutions. A new standalone version is being developed to address these constraints and provide more reliable notifications.

Thank you for your patience and understanding as we improve our service.

> **Lange Lab - Dec 2024**

## Introduction
The local sequence context is the most fundamental feature determining proteins' post-translational modification (PTM). Recent technological improvements allow for the detection of new and less prevalent modifications. We found that established state-of-the-art algorithms for the detection of PTM motifs in complex datasets failed to keep up with this technological development and are no longer robust. To overcome this limitation, we developed RoLiM, a new linear motif deconvolution algorithm and a web server that enables robust and unbiased identification of local amino acid sequence determinants in complex biological systems.


## How to use?
For convenient use, we provide a web frontend accessible at http://langelab.org/rolim

If preferred, RoLiM can also be installed locally, but the software has not been designed to ease local installation and use. 

The webfrontent includes detailed explanations for all options and parameters as well as example datasets for download.

`Example data 1:` prealigned sequence list (http://langelab.org/rolim/textfile)
```splus
HPKPKQFSSFEKRAK
DVATSPISPTENNTT
DLQEVLSSDENGGTY
EPDHYRYSDTTDSDP
TETRSSSSESSHSSS
GDDEDACSDTEATEA
DPEKFADSDQDRDPH
PEPSTKVSEEAESQQ
DVHMVSDSDGDDFED
EGASLELSDDDTESK
FLWSPFESLCEIGEK
NSYSGSNSGAAIGWG
GEYRSLESDNEEKAF
SSGSASKSDKDLETQ
EGRGEVGSAGDMRAA
EDLVDSLSEGDAYPN
EASHSGGSGDEAPKL
KGATKESSEKDRGRD
TPEELDDSDFETEDF
EVVGGDDSDGLRAED
KRSRAIHSSDEGEDQ
PASCPLDSDLSEDED
TDNLLPMSPEEFDEV
RSDVESSSEEEDVTT
DEEDLVDSLSEGDAY
SPDIDNYSEEEEESF
FKLFHPSSESEQGLT
SKGVDFESSEDDDDD
SMGLYMDSVRDADYS
STEQTLASDTDSSLD
PGETPPLSPIDMESQ
PGPQSPGSPLEEERQ
ESQDSSDSIGSSQKA
QKQEPLGSDSEGVNC
KDSSHYDSDGDKSDD
ASLSSLNSSESDKDQ
SADAANGSNEDRGEV
FHYRTLHSDDEGTVL
DDDDDDNSDEEDNDD
PGLQAADSDDEDLDD
IGDLVLDSDEEENGQ
IYPWMRSSGTDRKRG
DEELEGISPDELKDE
TSSYLSDSGSTGEHT
NSEASNASETESDHR
CPLDSDLSEDEDLQL
ITDVHMVSDSDGDDF
QMEKDIRSDVEESDS
PFAFNLNSDTDVEEG
QIRLRRDSKEANARR
```

`Example Data 2:` peptide list with optional protein identifiers (http://langelab.org/rolim/peptidelist)
```
VLYDVQELR	P09525
SLVINYDLPTNR	Q14240
DALGLNIYEQNDR	P26038
KLGIHEDSQNR	P07900
LATQLTGPVMPVR	P26373
ASSNESLVVNR	P42167
VVYGGADIGQQIR	O00571
AGDEIDEPSER	Q96ME7
SLGSALRPSTSR	P08670
QASQGTLQTR	P78527
NYGPMKSGNFGGSR	P22626
SELAGHQTSAESWGTGR	P36578
SSVSDFNQNFEVPNR	P41218
AAQSGILDR	Q9Y6N5
SGKASSAAGLTAAVVR	Q14566
GDIIIDGGNSEYR	P52209
SSLAEGSVTSVGSVNPAENFR	P13010
ADLQNDEVAFR	P61247
EVDNELR	Q9Y3Z3
GEKDIPGLTDTTVPR	P62753
EVDNELR	Q9Y3Z3
GEKDIPGLTDTTVPR	P62753
EGADNQGAGEQGR	P67809
MVQAEEAAAEITR	O00567
TEQGAELSNEER	P63104
QNSESGEKNEGSESAPEGQAQQR	P67809
SADDTPEVLNR	O95674
GGIKEDTEEHHLR	P09651
VSGVCVQTETVLR	P13639
NMACYCR	P59665
NMVPFPR	P07437
KAEPVEVVAPR	P09874
TDESLR	P09651
ALDTKGPEIR	P14618
SNLNPER	P26599
NDGEVDDEEDEEELGEEER	P39687
AASEFFR	P06733
SVVDLTCR	P04406
ILPPTRPTDKPLR	Q05639
VVCPKDEDYKQR	Q00839
VLYDVQELR	P09525
SLVINYDLPTNR	Q14240
DALGLNIYEQNDR	P26038
KLGIHEDSQNR	P07900
LATQLTGPVMPVR	P26373
ASSNESLVVNR	P42167
VVYGGADIGQQIR	O00571
AGDEIDEPSER	Q96ME7
SLGSALRPSTSR	P08670
QASQGTLQTR	P78527
```

- Download one of the provided demo datasets
- Enter your email, a short name and optional description and select the downloaded dataset as "Foreground dataset"
- Select the appropriate foreground dataset format
- Select the species for the background dataset ('Human')
- Leave all parameters in default setting and `Submit job`

## Results
Processing of the data by the server takes usually a few minutes or less and results as well as error messages will be sent to the provided email address. If the server is working on a long queue of submissions or the submitted dataset is particularly extensive or complex processing time can be considrably longer. 

Once completed you will receive an email with an .zip archive and the following content and folder structure:

```
- summary
-- log.txt  | summary of all parameters and options
-- sequence_clustermap.svg  | a hierarchically clustered heatmap of all provided sequences matched to the identified linear motifs / patterns
- patterns
-- pattern_summary_table.csv | tabular summary of all patterns and their respective sample and foregroundand frequencies and enrichement 
-- for each enriched pattern: 
-- pattern_sequences.txt | a list of all submitted sequences matching to a pattern
-- pattern_log_map.pdf | a weblogo representation of all submitted sequences matching to a pattern
```

## Citation
Currently under review

## Todo

- Add pre-computed backgrounds for species other than homo sapiens
- Add easy to interpret error messages in addition to the provided trace back. 
- Improve handling of some non-sensical parameter combinations to avoid failed submissions
- Improve support for Safari web browsers

## Feedback
If you have any feedback or suggestions, please reach out or open an issue on github.




## Dependencies:
```
amqp==2.4.1
backcall==0.1.0
billiard==3.5.0.5
celery==4.2.1
Click==7.0
cycler==0.10.0
decorator==4.3.2
Django==2.1.5
django-filter==2.1.0
django-rq==2.0
djangorestframework==3.9.1
gunicorn==19.9.0
ipykernel==5.1.0
ipython==7.2.0
ipython-genutils==0.2.0
jedi==0.13.2
jupyter-client==5.2.4
jupyter-core==4.4.0
kiwisolver==1.0.1
kombu==4.3.0
Markdown==3.0.1
matplotlib==3.0.2
matplotlib-venn==0.11.5
mysqlclient==1.4.2
numpy==1.16.1
pandas==0.24.1
parso==0.3.4
pexpect==4.6.0
pickleshare==0.7.5
prompt-toolkit==2.0.8
ptyprocess==0.6.0
Pygments==2.3.1
PyMySQL==0.9.3
pyparsing==2.3.1
python-dateutil==2.8.0
pytz==2018.9
pyzmq==17.1.2
redis==3.2.1
rq==1.0
scipy==1.2.0
seaborn==0.9.0
setuptools-scm==3.2.0
six==1.12.0
tornado==5.1.1
traitlets==4.3.2
vine==1.2.0
wcwidth==0.1.7
weblogo==3.7.1
```
