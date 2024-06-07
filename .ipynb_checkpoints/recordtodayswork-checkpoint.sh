#!/bin/bash

git pull
git add .
git add -u
now=$(date '+%d/%m/%Y %H:%M:%S')
git commit -m "$now"
git push
