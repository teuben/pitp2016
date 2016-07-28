#! /usr/bin/env bash

# requires pdftk (sudo apt-get install pdftk)

if [ -f AllLightningTalks.pdf ]; then
  rm AllLightningTalks.pdf
fi

pdftk *.pdf cat output AllLightningTalks.pdf
