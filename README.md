# smk-m2

Use snakemake to reproduce gatk best practice paired SNV calling [v2.4](https://github.com/gatk-workflows/gatk4-somatic-snvs-indels/tree/2.4.0). Due to TCGA legacy accession issue, we will first localize bams with signed URL.

*still under development!*

## Directory structure

```bash
├── cluster.json # configure cluster resources here
├── config.yaml # configure input params here
├── dag.svg
├── README.md
├── rules
│   └── localize.smk
├── scripts
│   ├── common
│   │   └── __init__.py
│   └── localize_from_nci.py
├── Snakefile
└── test_samples.txt # two tumor-normal pairs for testing
```

Please change the `credential` path from `config.yaml`to your own =)

## DAG

(splitted intervals into 3 segments for clearer viz)

![dag](https://user-images.githubusercontent.com/30106174/68724663-230e8000-058a-11ea-9104-ad8409e8588f.png)