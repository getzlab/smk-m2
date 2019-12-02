#!/bin/bash

sudo munged -f
sudo slurmd -f /demo-mount/slurm.conf
cp /demo-mount/slurm.conf /opt/etc/slurm.conf


echo 'export PATH=/demo-mount/gatk-4.1.4.0:$PATH' >> /home/qing/.bashrc && source /home/qing/.bashrc
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash Miniconda3-latest-Linux-x86_64.sh -b 
bash Miniconda3-latest-Linux-x86_64.sh -p /home/qing/miniconda3 -b

echo 'export PATH=/home/qing/miniconda3/bin:/home/qing/miniconda3/.local/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
conda install -y matplotlib numpy graphviz pandas && pip install snakemake 
chown qing -R /home/qing/miniconda3/

ln -s /home/qing/miniconda3/bin/python3 /usr/local/bin/python3
slurmd -f /demo-mount/slurm.conf
