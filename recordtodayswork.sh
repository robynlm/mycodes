#!/bin/bash

git pull
git add --all .
git add -u
now=$(date)
git commit -m "$now"
git push
