#!/bin/bash

function helpTxt(){
	echo -e "\nDESCRIPTION:\n  Wrapper script for running ICA_AROMA through our singularity container. See the ICA_AROMA documentation at https://github.com/maartenmennes/ICA-AROMA for more robust descriptions\nUSAGE:\n  ./$(basename $0) [OPTIONS]\nOPTIONS:\n  -i       Input nifti file path\n  -o       Output directory\n  -mc      realignment parameters text file\n  -tr      repetition time of functional image (in seconds)\n  -den     denoising strategy\n  -m       mask file for MELODIC\n  -h       Display this help text\n"
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
		in="$2"
		shift 2
		;;
    -o)
    out="$2"
    shift 2
    ;;
    -mc)
    mc="$2"
    shift 2
    ;;
    -tr)
    tr="$2"
    shift 2
    ;;
    -den)
    den="$2"
    shift 2
    ;;
    -m)
    mask="$2"
    shift 2
    ;;
    -h)
    helpTxt
    ;;
	esac
done

export FSLOUTPUTTYPE="NIFTI_GZ"

/opt/ICA-AROMA/ICA_AROMA.py -i ${in} -o ${out} -mc ${mc} -tr ${tr} -den ${den} -m ${mask}
