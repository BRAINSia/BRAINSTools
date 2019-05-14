#!/bin/bash

$(dirname "$0")/DevelopmentSetupScripts/SetupHooks.sh
## Add ssh commit reference
git remote add bia git@github.com:BRAINSia/BRAINSTools.git
## run pre-commit for python dir AutoWorkup
pip install -r $(dirname "$0")/../AutoWorkup/requirements-dev.txt
pre-commit install
