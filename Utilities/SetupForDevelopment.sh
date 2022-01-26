#!/bin/bash

$(dirname "$0")/DevelopmentSetupScripts/SetupHooks.sh
## Add ssh commit reference
git remote -v | grep bia || git remote add bia git@github.com:BRAINSia/BRAINSTools.git 2>&1 || true
git fetch bia
## run pre-commit for python dir AutoWorkup
#pip install -r $(dirname "$0")/../AutoWorkup/requirements-dev.txt
#pre-commit install
