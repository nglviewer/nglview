#!/bin/sh

if [[ "$PYTHON_VERSION" == "3.5" ]]; then
    sleep 3
    jupyter notebook --port=8889 &
    sleep 3
    nightwatch
fi
