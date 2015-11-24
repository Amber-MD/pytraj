#!/bin/sh

if [[ "$PYTHON_VERSION" != "3.4" ]]; then
    coveralls
fi
