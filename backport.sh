#!/bin/bash
# Script to strip f-strings out of frbpoppy and replace with the older .format
# Necessary for older Python installations of <=3.5

pip3 install future-fstrings[rewrite]
for e in ./frbpoppy ./tests; do
  cd $e
  for f in *.py; do
    echo $f
    tmpfile=$(mktemp)
    future-fstrings-show $f > ${tmpfile}
    cat ${tmpfile} > $f
    rm -f ${tmpfile}
  done
  cd ..
done
