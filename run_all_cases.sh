#!/bin/bash 

set -euxo pipefail

python run_simulation.py -k 10000 -D 100 -e 0.005 -f 0.01 -a 0.05 -t 2 -s 42
python run_simulation.py -k 10000 -D 300 -e 0.005 -f 0.01 -a 0.05 -t 2 -s 42
python create_af_plots_for_simulation.py -k 10000 -D 300 -e 0.005 -f 0.01 -a 0.2 0.1 0.05 0.02 0.01 0.005 -t 2 -s 42
python create_dp_plots_for_simulation.py -k 10000 -D 10 20 30 40 50 100 150 200 210 220 230 240 250 260 270 280 290 300 310 320 330 340 350 360 370 380 390 400 -e 0.005 -f 0.01 -a 0.05 -t 2 -s 42
python run_simulation.py -k 3000000000 -D 300 -e 0.005 -f 0.01 -a 0.05 -t 2 -s 42

echo "Done!"
