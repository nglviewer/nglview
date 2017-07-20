#!/bin/sh

echo "Running notebook in background"
jupyter notebook --port=8889 > /dev/null &
sleep 2
