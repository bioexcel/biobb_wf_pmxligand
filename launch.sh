#!/bin/bash

source env.sh

enqueue_compss --num_nodes=2 --exec_time=120 --network=ethernet --worker_working_dir=$PWD --master_working_dir=$PWD --base_log_dir=$PWD \
 --qos=debug -d -g --worker_in_master_cpus=48 --jvm_workers_opts="-Dcompss.worker.removeWD=false" \
 /home/bsc23/bsc23210/macshare/biobb_wf_pmxligand/pmxLig.py --config /home/bsc23/bsc23210/macshare/biobb_wf_pmxligand/pmxLig.yaml

