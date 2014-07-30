#!/bin/bash

$(dirname "$0")/DevelopmentSetupScripts/SetupHooks.sh
## Add ssh commit reference
git remote add bia git@github.com:BRAINSia/BRAINSTools.git
