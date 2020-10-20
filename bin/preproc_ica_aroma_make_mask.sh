#!/bin/bash

function helpTxt(){
	echo -e "\nDESCRIPTION:\n  Make a brain extracted mask of mean functional image for use in ICA AROMA with SPM data.\nUSAGE:\n  ./$(basename $0) -i <inputNii> -o <outputMask>\nOPTIONS:\n  -i       Input nifti file path\n  -o       Output mask file path\n  -h       Display this help text\n"
	exit
}

#========================================
#========parse/validate inputs===========
#========================================
while [[ $# -gt 0 ]]
do
	key="$1"
	case $key in
		-i)
		niiFile="$2"
		shift 2
		;;
    -o)
    mask="$2"
    shift 2
    ;;
    -h)
    helpTxt
    ;;
	esac
done

export FSLOUTPUTTYPE="NIFTI_GZ"

#make tmp filepath for mean func image
meanFunc="$(dirname "${niiFile}")/mean_$(basename "${niiFile}")"

#make mean func image
fslmaths "${niiFile}" -Tmean "${meanFunc}"

#run bet on mean func image. FSL takes the output filename and appends a "_mask"
#before the extension, so we'll need to move the file to the requested filename
#after running
bet "${meanFunc}" "${mask}" -f 0.3 -n -m -R
betOut="${mask/".nii.gz"/}_mask.nii.gz"
mv "${betOut}" "${mask}"
