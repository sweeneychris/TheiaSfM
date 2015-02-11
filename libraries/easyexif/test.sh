#!/bin/bash
for jpeg in `ls test-images/*.jpg`; do
  ./demo $jpeg > /tmp/`basename $jpeg`.actual
  diff $jpeg.expected /tmp/`basename $jpeg`.actual > /tmp/diff.out
  if [[ -s /tmp/diff.out ]] ; then
    echo "FAILED ON $jpeg"
    cat /tmp/diff.out
    exit
  fi ;
  echo "PASS $jpeg"

done
