#!/bin/sh

echo "Checking ./nbtests"

if [ ! -d nbtests ]; then
    echo "Can not find nbtests folder"
    echo "Clone it"
    git clone https://github.com/hainm/nbtests
fi

echo "Running notebook in background"
jupyter notebook --port=8889 --NotebookApp.iopub_data_rate_limit 1000000000 > /dev/null &
sleep 2
