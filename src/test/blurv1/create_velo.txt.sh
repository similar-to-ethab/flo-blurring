#!/bin/bash

cd $BLUR/blur

ls -v *-velo-*.png | awk '{print "file \x27png-"$0"\x27"}' >> velo.txt

