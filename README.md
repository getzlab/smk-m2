# smk-m2

Use snakemake to reproduce gatk best practice paired SNV calling [v2.4](https://github.com/gatk-workflows/gatk4-somatic-snvs-indels/tree/2.4.0). Due to TCGA legacy accession issue, we will first localize bams with signed URL.

*still under development!*

## Directory structure

```bash
.
├── cluster.json
├── config.yaml
├── dag.svg
├── README.md
├── res
├── scripts
│   ├── common
│   │   └── __init__.py
│   └── localize_from_nci.py
├── Snakefile
└── test_samples.txt
```

Please change the `credential` path from `config.yaml`to your own =)

## DAG

(splitted intervals into 3 segments for clearer viz)

![dag](dag.svg)