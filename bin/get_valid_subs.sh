#!/bin/bash

#read in entries from the sm_config file
task=$(cat config/sm_config.json | jq -r ".task[]")
ses=$(cat config/sm_config.json | jq -r ".ses[]")
acq=$(cat config/sm_config.json | jq -r ".acq[]")
run=$(cat config/sm_config.json | jq -r ".run[]")
bidsDir=$(cat config/sm_config.json | jq -r ".bids_dir")

#expected input files
reqd='${bidsDir}/sub-${sub}/ses-${s}/func/sub-${sub}_ses-${s}_task-${t}_acq-${acq}_run-${r}_bold.nii.gz ../../reorientation/sub-${sub}/ses-${s}/T1w_reorient.mat ${bidsDir}/sub-${sub}/ses-${s}/anat/sub-${sub}_ses-${s}_acq-highres_T1w.nii.gz'

#initialize empty subs array
subs=""
echo "${bidsDir}"
for sub in "${bidsDir}"/sub-*; do
  #extract subject id from filepath
  sub=$(echo ${sub} | cut -d "/" -f 4 | cut -d "-" -f 2)
  addSub=1
  missing="subject ${sub} is missing the following required files:"
  for t in ${task}; do
    for s in ${ses}; do
      for r in ${run}; do
        for f in ${reqd}; do
          f=$(eval echo "${f}")
          if ! [ -f "${f}" ]; then
            addSub=0
            missing=${missing}$'\n'${f}
          fi
        done
      done
    done
  done
  if [ ${addSub} -eq 1 ]; then
    subs="${subs} \"${sub}\","
  else
    echo "${missing}"
    echo
  fi
done

echo
echo "${subs}"
