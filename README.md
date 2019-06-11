# CNV caller

## Installation
```
python3 -m venv venv
. ./venv/bin/activate
pip install -r requirements.txt
```

## Usage
```
python cnv_caller/call.py --verbose --tumour tumour.bam --normal normal.bam --bed regions.bed > cnv.tsv # slow step
python cnv_caller/group.py --minbases 1000 --minlen 10000 < cnv.tsv > cnvs.tsv # fast step
```

## Method
* call.py measures the ratio of depths between regions specified by the bed file
* group.py decides if a group of regions form a CNV, using the minbases and minlen parameters
