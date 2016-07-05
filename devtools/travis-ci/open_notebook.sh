#!/bin/sh

if [[ "$PYTHON_VERSION" == "3.5" ]]; then
    sleep 3
    sudo jupyter notebook --port=8889 --browser=google-chrome &
    sleep 3
    nightwatch
fi
