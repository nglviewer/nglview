#!/bin/sh

echo "Running notebook in background"
jupyter notebook --port=8889 --NotebookApp.iopub_data_rate_limit 1000000000 > /dev/null &
sleep 2
