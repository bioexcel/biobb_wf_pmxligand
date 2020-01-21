#!/bin/bash

source env.sh

enqueue_compss --num_nodes=1 --exec_time=10 --network=ethernet --worker_working_dir=$PWD \
 --qos=debug -d --worker_in_master_cpus=48 --jvm_workers_opts="-Dcompss.worker.removeWD=false" \
 pmxLig.py --config pmxLig.yaml

