#!/bin/bash

cd build
rm -fr *
mkdir bin
cp ../SetUpData.ini bin/
cp -r ../Mesh/ bin/
cd ../