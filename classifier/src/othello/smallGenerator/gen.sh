#!/bin/bash
./gen 1500 20 128 | sort > smallkmers.txt
./gen 15000000 20 128 | sort > 15M.txt
