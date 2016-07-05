#!/bin/sh

if [[ "$PYTHON_VERSION" == "3.5" && "$TRAVIS_OS_NAME" == "osx" ]]; then
    sleep 3
    jupyter notebook --port=8889 --browser=google-chrome &
    sleep 3
    nightwatch
fi
