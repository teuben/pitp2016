#! /usr/bin/env bash

if [ -f AllLightningTalks.pdf ]; then
  rm AllLightningTalks.pdf
fi

pdftk *.pdf cat output AllLightningTalks.pdf
