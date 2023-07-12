#!/bin/bash

env_name=somatic_snv_sim
conda create -n $env_name python=3.10
conda activate $env_name
pip install -r requirements.txt
