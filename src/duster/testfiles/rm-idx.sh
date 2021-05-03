#! /bin/bash


for i in *.kidx; do
  if [ -f $i ]; then
    rm $i
  fi
done