#!/bin/sh

platform=`python -c "import sys; print(sys.platform)"`
if [ $platform != 'linux' ]; then
    echo "do nothing, let pytraj handle"
else
   python scripts/install_libcpptraj.py
fi
