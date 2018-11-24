#!/bin/bash
# \author Hans J. Johnson
#
# Push data to new repo with sha512 listing


for md5file in *.md5; do
  sha512file=${md5file//.md5/.sha512}
  realfile=${md5file//.md5/}
  if [[ -f ${realfile} ]]; then
     ~/Dashboard/src/ITK/Utilities/UploadBinaryData.sh --folder-id 5be8b0308d777f2179a18270 ${realfile}
     if [[ $? -eq 0 ]]; then
        git add -- ${sha512file}
        git rm  -- ${md5file}
     fi
  fi
done
