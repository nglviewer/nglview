#!/bin/sh

if [[ "$PYTHON_VERSION" == "3.5" ]]; then
    sleep 3
    sudo /home/travis/miniconda/envs/myenv/bin/jupyter notebook --port=8889 --browser=google-chrome --ip=127.0.0.1 &
    sleep 3
    nightwatch
fi
