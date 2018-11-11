#!/bin/bash
# \author Hans J. Johnson

old_reference=$1
new_reference=$2

old_basename=$(basename ${old_reference})
old_hash_file=$(find ./TestData -name ${old_basename}.*)

if [[ -f ${old_hash_file} ]]; then
  git rm ${old_hash_file}
fi
if [[ -f ${old_reference} ]]; then
  rm ${old_reference}
fi

new_hash_file_temp=$(echo ${old_hash_file} | sed 's/.md5//g' |sed 's/.sha512//g' )

cp ${new_reference} ${new_hash_file_temp}
~/Dashboard/src/ITK/Utilities/UploadBinaryData.sh --folder-id 5be8b0308d777f2179a18270 ${new_hash_file_temp}
git add -- ${new_hash_file_temp}.sha512
