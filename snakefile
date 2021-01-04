import os, subprocess

localrules: all, getRawFunc, getRawT1w, gunzipIcaAroma, covAAL90, covAAL116, covDiedrichsen2011, covGordon2016, covPower264, covSchaefer2018, covTian2020Subcortical, symBold, symArt, symRp, symT1w
configfile: 'config/sm_config.json'

rule all:
    input:
        expand('data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/art_regression_outliers_swrtrun-{run}.mat',
             sub=config['sub'],ses=config['ses'],acq=config['acq'],task=config['task'],run=config['run']),
        expand('data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-{run}_icaaroma/denoised_func_data_nonaggr.nii',
            sub=config['sub'],ses=config['ses'],acq=config['acq'],task=config['task'],run=config['run']),
        expand('data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-{run}/cov/AAL90',
            sub=config['sub'],ses=config['ses'],acq=config['acq'],task=config['task'],run=config['run']),
        expand('data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-{run}/cov/AAL116',
            sub=config['sub'],ses=config['ses'],acq=config['acq'],task=config['task'],run=config['run']),
        expand('data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-{run}/cov/Diedrichsen2011',
            sub=config['sub'],ses=config['ses'],acq=config['acq'],task=config['task'],run=config['run']),
        expand('data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-{run}/cov/Gordon2016',
            sub=config['sub'],ses=config['ses'],acq=config['acq'],task=config['task'],run=config['run']),
        expand('data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-{run}/cov/Power264',
            sub=config['sub'],ses=config['ses'],acq=config['acq'],task=config['task'],run=config['run']),
        expand('data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-{run}/cov/Schaefer2018',
            sub=config['sub'],ses=config['ses'],acq=config['acq'],task=config['task'],run=config['run']),
        expand('data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-{run}/cov/Tian2020',
            sub=config['sub'],ses=config['ses'],acq=config['acq'],task=config['task'],run=config['run']),
        expand('../pipeline-resting-L1/data/preprocessed/sub-{sub}/ses-{ses}/sub-{sub}_ses-{ses}_task-{task}_acq-{acq}_run-{run}_bold.nii',
            sub=config['sub'],ses=config['ses'],acq=config['acq'],task=config['task'],run=config['run']),
        expand('../pipeline-resting-L1/data/preprocessed/sub-{sub}/ses-{ses}/sub-{sub}_ses-{ses}_space-{task}.{acq}_run-{run}_hT1w.nii',
            sub=config['sub'],ses=config['ses'],acq=config['acq'],task=config['task'],run=config['run']),
        expand('../pipeline-resting-L1/data/preprocessed/sub-{sub}/ses-{ses}/sub-{sub}_ses-{ses}_task-{task}_acq-{acq}_run-{run}_art_regression_outliers.mat',
            sub=config['sub'],ses=config['ses'],acq=config['acq'],task=config['task'],run=config['run']),
        expand('../pipeline-resting-L1/data/preprocessed/sub-{sub}/ses-{ses}/sub-{sub}_ses-{ses}_task-{task}_acq-{acq}_run-{run}_rp.txt',
            sub=config['sub'],ses=config['ses'],acq=config['acq'],task=config['task'],run=config['run'])

# STEP 1 - COPY INPUT FILES
rule getRawFunc:
    input:
        config['bids_dir']+'/sub-{sub}/ses-{ses}/func/sub-{sub}_ses-{ses}_task-{task}_acq-{acq}_run-{run}_bold.nii.gz'
    output:
        'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-{run}.nii'
    shell:
        'gunzip -c {input} > {output}'

rule getRawT1w:
    input:
        RawT1w = config['bids_dir']+'/sub-{sub}/ses-{ses}/anat/sub-{sub}_ses-{ses}_acq-highres_T1w.nii.gz'
    output:
        t1w = temp('data/anat/sub-{sub}/ses-{ses}/T1w.nii')
    shell:
        'gunzip -c {input} > {output}'

# STEP 2 - PREPROCESS

rule hCorr:
	input:
	   t1w = 'data/anat/sub-{sub}/ses-{ses}/T1w.nii'
	output:
		hT1w = 'data/anat/sub-{sub}/ses-{ses}/hT1w.nii'
	shell:
		'''
		matlab -nodisplay -r "cd $PWD; addpath $PWD/bin/; addpath {config[spm_dir]};
		hCorr('{input.t1w}','{output.hT1w}'); exit;"
		'''.replace('\n','')

rule Unifize:
    input:
        'data/anat/sub-{sub}/ses-{ses}/hT1w.nii'
    output:
        'data/anat/sub-{sub}/ses-{ses}/uhT1w.nii'
    shell:
        '''
        singularity run -B $PWD:/data ../../lib/fmriproc.sif 3dUnifize -prefix /data/{output} /data/{input}
        '''

#For spiral only
#append functional image so 10th volume of ref image is first (and thus used as reference)
rule prepRealign:
    input:
        refFunc = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-1.nii',
        func = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-{run}.nii'
    output:
        'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/TMP_run-{run}.nii'
    shell:
        '''
        if [ -e {input.func}.gz ]; then rm {input.func}.gz; fi
        FSLOUTPUTTYPE=NIFTI
        funcPath=$(dirname "{input.refFunc}")
        refVol=$funcPath/tmpRefVol_{wildcards.run}
        echo "extracting 10th volume from reference functional..."
        fslroi {input.refFunc} $refVol 9 1 #align to raw 10th volume
        tmpStr=$funcPath/tmp_tFunc_{wildcards.run}_vol_
        echo "appending 10th volume to temporary functional nifti..."
        fslsplit {input.func} $tmpStr -t
        fslmerge -t {output} $refVol $tmpStr*
        rm $refVol* $tmpStr*
        '''

rule spmRealign:
    input:
        'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/TMP_run-{run}.nii'
    output:
        tmpAFunc = temp('data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/rTMP_run-{run}.nii'),
        tmpRpTxt = temp('data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/rp_TMP_run-{run}.txt')
    shell:
        '''
        tmpJob=$(mktemp tmp/XXXXXXXX.mat)
        matlab -nodisplay -r "cd $PWD; addpath $PWD/bin; addpath {config[spm_dir]}; \
        preproc_spm_realign('inFunc', '{input}','outFile', '$tmpJob'); exit"
        singularity run -B $PWD:/data {config[spmSif]} /data/bin/spm_jobman_run.m $tmpJob
        rm $tmpJob
        '''

#remove 1st volume of realign output (previously appended)
#remove 1st entry to rpTxt file, which corresponded with ref image
rule cleanRealign:
    input:
        tmpRFunc = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/rTMP_run-{run}.nii',
        tmpRpTxt = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/rp_TMP_run-{run}.txt'
    output:
        rFunc = temp('data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/rrun-{run}.nii'),
        rpTxt = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/rp_run-{run}.txt'
    shell:
        '''
        FSLOUTPUTTYPE=NIFTI
        fslroi {input.tmpRFunc} {output.rFunc} 1 -1
        sed -e "1d" {input.tmpRpTxt} > {output.rpTxt}
        '''

rule spmSliceTime:
    input:
        func = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-{run}.nii',
        json = config['bids_dir']+'/task-{task}_acq-{acq}_bold.json'
    output:
        trFunc = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/trun-{run}.nii',
    shell:
        '''
        tmpJob=$(mktemp tmp/XXXXXXXX.mat)
        matlab -nodisplay -r "cd $PWD; addpath $PWD/bin; addpath {config[spm_dir]}; \
        preproc_spm_slicetime('inFunc','{input.func}','outFile','$tmpJob', \
        'bidsConfig', '{input.json}'); exit"
        singularity run -B $PWD:/data {config[spmSif]} /data/bin/spm_jobman_run.m $tmpJob
        rm $tmpJob
        '''

# now we realign again, this time with the slicetime corrected functional
rule prepRealign2:
    input:
        refFunc = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-1.nii',
        tFunc = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/trun-{run}.nii'
    output:
        'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/TMP_trun-{run}.nii'
    shell:
        '''
        if [ -e {input.tFunc}.gz ]; then rm {input.tFunc}.gz; fi
        FSLOUTPUTTYPE=NIFTI
        funcPath=$(dirname "{input.refFunc}")
        refVol=$funcPath/tmpRefVol_{wildcards.run}
        echo "extracting 10th volume from reference functional..."
        fslroi {input.refFunc} $refVol 9 1 #align to raw 10th volume
        tmpStr=$funcPath/tmp_tFunc_{wildcards.run}_vol_
        echo "appending 10th volume to temporary functional nifti..."
        fslsplit {input.tFunc} $tmpStr -t
        fslmerge -t {output} $refVol $tmpStr*
        rm $refVol* $tmpStr*
        '''

rule spmRealign2:
    input:
        'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/TMP_trun-{run}.nii'
    output:
        tmpAFunc = temp('data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/rTMP_trun-{run}.nii'),
        tmpRpTxt = temp('data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/rp_TMP_trun-{run}.txt')
    shell:
        '''
        tmpJob=$(mktemp tmp/XXXXXXXX.mat)
        matlab -nodisplay -r "cd $PWD; addpath $PWD/bin; addpath {config[spm_dir]}; \
        preproc_spm_realign('inFunc', '{input}','outFile', '$tmpJob'); exit"
        singularity run -B $PWD:/data {config[spmSif]} /data/bin/spm_jobman_run.m $tmpJob
        rm $tmpJob
        '''

#remove 1st volume of realign output (previously appended)
#remove 1st entry to rpTxt file, which corresponded with ref image
rule cleanRealign2:
    input:
        tmpRFunc = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/rTMP_trun-{run}.nii',
        tmpRpTxt = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/rp_TMP_trun-{run}.txt'
    output:
        rtFunc = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/rtrun-{run}.nii',
        rpTxt = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/rp_trun-{run}.txt'
    shell:
        '''
        FSLOUTPUTTYPE=NIFTI
        fslroi {input.tmpRFunc} {output.rtFunc} 1 -1
        sed -e "1d" {input.tmpRpTxt} > {output.rpTxt}
        '''

rule spmCoregister:
    input:
        rtFunc = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/rtrun-1.nii',
        uhT1w = 'data/anat/sub-{sub}/ses-{ses}/uhT1w.nii'
    output:
        cuhT1w = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/cuhT1w.nii'
    shell:
        '''
        tmpJob=$(mktemp tmp/XXXXXXXX.mat)
        cp {input.uhT1w} {output.cuhT1w}
        preTime=$(stat -c %y {output.cuhT1w})
        matlab -nodisplay -r "cd $PWD; addpath $PWD/bin; addpath {config[spm_dir]}; \
        preproc_spm_coregister('refImg', '{input.rtFunc},10', 'srcImg', '{output.cuhT1w}', \
        'outFile', '$tmpJob'); exit"
        singularity run -B $PWD:/data {config[spmSif]} /data/bin/spm_jobman_run.m $tmpJob
        postTime=$(stat -c %y {output.cuhT1w})
        if [[ "$preTime" == "$postTime" ]]; then
            rm {output.cuhT1w}
        fi
        rm $tmpJob
        '''

rule spmSegment:
    input:
        cuhT1w = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/cuhT1w.nii',
        TPM = 'lib/tpm/SPM_TPM.nii'
    output:
        defField = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/y_cuhT1w.nii',
        gm = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/wc1cuhT1w.nii',
        wm = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/wc2cuhT1w.nii',
        csf = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/wc3cuhT1w.nii'
    shell:
        '''
        tmpJob=$(mktemp tmp/XXXXXXXX.mat)
        matlab -nodisplay -r "cd $PWD; addpath $PWD/bin; addpath {config[spm_dir]}; \
        preproc_spm_segment('inStruct','{input.cuhT1w}', \
        'TPM', '{input.TPM}', 'outFile', '$tmpJob'); exit"
        singularity run -B $PWD:/data {config[spmSif]} /data/bin/spm_jobman_run.m $tmpJob
        rm $tmpJob
        '''

rule spmNormalizeFunc:
    input:
        rtFunc = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/rtrun-{run}.nii',
        defField = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/y_cuhT1w.nii'
    output:
        wrtFunc = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/wrtrun-{run}.nii',
    shell:
        '''
        tmpJob=$(mktemp tmp/XXXXXXXX.mat)
        matlab -nodisplay -r "cd $PWD; addpath $PWD/bin; addpath {config[spm_dir]}; \
        preproc_spm_normalize_write('inDef', '{input.defField}', 'inNii', '{input.rtFunc}', \
        'outFile', '$tmpJob'); exit"
        singularity run -B $PWD:/data {config[spmSif]} /data/bin/spm_jobman_run.m $tmpJob
        rm $tmpJob
        '''

rule spmNormalizeT1w:
    input:
        cuhT1w = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/cuhT1w.nii',
        defField = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/y_cuhT1w.nii'
    output:
        wcuhT1w = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/wcuhT1w.nii',
    shell:
        '''
        tmpJob=$(mktemp tmp/XXXXXXXX.mat)
        matlab -nodisplay -r "cd $PWD; addpath $PWD/bin; addpath {config[spm_dir]}; \
        preproc_spm_normalize_write('inDef', '{input.defField}', 'inNii', '{input.cuhT1w}', \
        'outFile', '$tmpJob'); exit"
        singularity run -B $PWD:/data {config[spmSif]} /data/bin/spm_jobman_run.m $tmpJob
        rm $tmpJob
        '''

rule spmSmooth:
    input:
        wrtFunc = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/wrtrun-{run}.nii'
    output:
        swrtFunc = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/swrtrun-{run}.nii'
    shell:
        '''
        tmpJob=$(mktemp tmp/XXXXXXXX.mat)
        matlab -nodisplay -r "cd $PWD; addpath $PWD/bin; addpath {config[spm_dir]}; \
        preproc_spm_smooth('inFunc', '{input.wrtFunc}', 'outFile', '$tmpJob'); exit"
        singularity run -B $PWD:/data {config[spmSif]} /data/bin/spm_jobman_run.m $tmpJob
        rm $tmpJob
        '''

rule mkArtConfig:
    input:
        swutFunc = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/swrtrun-{run}.nii',
        rpTxt = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/rp_run-{run}.txt'
    output:
        artCfg = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/art_run-{run}.cfg'
    script:
        './bin/preproc_art_mkconfig.py'

rule ART:
    input:
        artCfg = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/art_run-{run}.cfg'
    output:
        art = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/art_regression_outliers_swrtrun-{run}.mat'
    shell:
        '''
        matlab -nodisplay -r "cd $PWD; addpath $PWD/bin/; addpath
        $PWD/lib/art; addpath {config[spm_dir]};
        preproc_art_run('{input.artCfg}'); exit;"
        '''.replace('\n','')

rule icaAROMA:
    input:
        swrtFunc = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/swrtrun-{run}.nii',
        rpTxt = 'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/rp_run-{run}.txt',
        bidsJson = config['bids_dir']+'/task-{task}_acq-{acq}_bold.json',
        mask = 'lib/masks/MNI152_T1_2mm_mask.nii.gz'
    output:
        temp('data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-{run}_icaaroma/denoised_func_data_nonaggr.nii.gz')
    threads: 2
    shell:
        '''
        tr=$(cat {input.bidsJson} | jq -r '.RepetitionTime')
        singularity run -B $PWD:/data ../../lib/fmriproc.sif /data/bin/preproc_ica_aroma_run.sh                       \
            -i /data/{input.swrtFunc} -o /data/$(dirname "{output}") -mc /data/{input.rpTxt}    \
            -tr $tr -den both -m /data/{input.mask} || true
        cd $(dirname {output})
        '''

rule unzipIcaAROMA:
    input:
        'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-{run}_icaaroma/denoised_func_data_nonaggr.nii.gz'
    output:
        'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-{run}_icaaroma/denoised_func_data_nonaggr.nii'
    shell:
        'gunzip {input}'


### COVERAGE CHECKING

def get_masks(masks):
    return [ os.path.join(masks,x) for x in os.listdir(masks) if '.nii' in x ]

rule covAAL90:
    input:
        nii='data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-{run}_icaaroma/denoised_func_data_nonaggr.nii',
        masks=get_masks(config['masks']['AAL90'])
    output:
        covDir=directory('data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-{run}/cov/AAL90')
    shell:
        '''
        matlab -nodisplay -r "cd $PWD; addpath $PWD/bin/; addpath
        $PWD/{config[spm_dir]}; batch_check_coverage('{input.nii}',
        '{input.masks}','{output.covDir}'); exit;"
        '''.replace('\n','')

rule covAAL116:
    input:
        nii='data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-{run}_icaaroma/denoised_func_data_nonaggr.nii',
        masks=get_masks(config['masks']['AAL116'])
    output:
        covDir=directory('data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-{run}/cov/AAL116')
    shell:
        '''
        matlab -nodisplay -r "cd $PWD; addpath $PWD/bin/; addpath
        $PWD/{config[spm_dir]}; batch_check_coverage('{input.nii}',
        '{input.masks}','{output.covDir}'); exit;"
        '''.replace('\n','')

rule covDiedrichsen2011:
    input:
        nii='data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-{run}_icaaroma/denoised_func_data_nonaggr.nii',
        masks=get_masks(config['masks']['Diedrichsen2011'])
    output:
        covDir=directory('data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-{run}/cov/Diedrichsen2011')
    shell:
        '''
        matlab -nodisplay -r "cd $PWD; addpath $PWD/bin/; addpath
        $PWD/{config[spm_dir]}; batch_check_coverage('{input.nii}',
        '{input.masks}','{output.covDir}'); exit;"
        '''.replace('\n','')

rule covGordon2016:
    input:
        nii='data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-{run}_icaaroma/denoised_func_data_nonaggr.nii',
        masks=get_masks(config['masks']['Gordon2016'])
    output:
        covDir=directory('data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-{run}/cov/Gordon2016')
    shell:
        '''
        matlab -nodisplay -r "cd $PWD; addpath $PWD/bin/; addpath
        $PWD/{config[spm_dir]}; batch_check_coverage('{input.nii}',
        '{input.masks}','{output.covDir}'); exit;"
        '''.replace('\n','')

rule covPower264:
    input:
        nii='data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-{run}_icaaroma/denoised_func_data_nonaggr.nii',
        masks=get_masks(config['masks']['Power264'])
    output:
        covDir=directory('data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-{run}/cov/Power264')
    shell:
        '''
        matlab -nodisplay -r "cd $PWD; addpath $PWD/bin/; addpath
        $PWD/{config[spm_dir]}; batch_check_coverage('{input.nii}',
        '{input.masks}','{output.covDir}'); exit;"
        '''.replace('\n','')

rule covSchaefer2018:
    input:
        nii='data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-{run}_icaaroma/denoised_func_data_nonaggr.nii',
        masks=get_masks(config['masks']['Schaefer2018'])
    output:
        covDir=directory('data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-{run}/cov/Schaefer2018')
    shell:
        '''
        matlab -nodisplay -r "cd $PWD; addpath $PWD/bin/; addpath
        $PWD/{config[spm_dir]}; batch_check_coverage('{input.nii}',
        '{input.masks}','{output.covDir}'); exit;"
        '''.replace('\n','')

rule covTian2020Subcortical:
    input:
        nii='data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-{run}_icaaroma/denoised_func_data_nonaggr.nii',
        masks=get_masks(config['masks']['Tian2020'])
    output:
        covDir=directory('data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-{run}/cov/Tian2020')
    shell:
        '''
        matlab -nodisplay -r "cd $PWD; addpath $PWD/bin/; addpath
        $PWD/{config[spm_dir]}; batch_check_coverage('{input.nii}',
        '{input.masks}','{output.covDir}'); exit;"
        '''.replace('\n','')

### MAKE SYMLINKS
rule symBold:
    input:
        'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/run-{run}_icaaroma/denoised_func_data_nonaggr.nii'
    output:
        '../pipeline-resting-L1/data/preprocessed/sub-{sub}/ses-{ses}/sub-{sub}_ses-{ses}_task-{task}_acq-{acq}_run-{run}_bold.nii'
    priority: 100
    shell:
        '''
        mkdir -p $(dirname {output})
        linkPath=$(realpath --relative-to=$(dirname {output}) $(dirname {input}))
        ln -s $linkPath/$(basename {input}) {output}
        '''

rule symT1w:
    input:
        'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/chT1w.nii'
    output:
        '../pipeline-resting-L1/data/preprocessed/sub-{sub}/ses-{ses}/sub-{sub}_ses-{ses}_space-{task}.{acq}_run-{run}_hT1w.nii'
    shell:
        '''
        mkdir -p $(dirname {output})
        linkPath=$(realpath --relative-to=$(dirname {output}) $(dirname {input}))
        ln -s $linkPath/$(basename {input}) {output}
        '''

rule symArt:
    input:
        'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/art_regression_outliers_swrtrun-{run}.mat'
    output:
        '../pipeline-resting-L1/data/preprocessed/sub-{sub}/ses-{ses}/sub-{sub}_ses-{ses}_task-{task}_acq-{acq}_run-{run}_art_regression_outliers.mat'
    shell:
        '''
        mkdir -p $(dirname {output})
        linkPath=$(realpath --relative-to=$(dirname {output}) $(dirname {input}))
        ln -s $linkPath/$(basename {input}) {output}
        '''


rule symRp:
    input:
        'data/acq-{acq}/preproc/task-{task}/sub-{sub}/ses-{ses}/rp_run-{run}.txt'
    output:
        '../pipeline-resting-L1/data/preprocessed/sub-{sub}/ses-{ses}/sub-{sub}_ses-{ses}_task-{task}_acq-{acq}_run-{run}_rp.txt'
    shell:
        '''
        mkdir -p $(dirname {output})
        linkPath=$(realpath --relative-to=$(dirname {output}) $(dirname {input}))
        ln -s $linkPath/$(basename {input}) {output}
        '''


'''
# 80% non-spiked volumes from ART
# Coverage Checks: Each ROI
1) 50% in each ROI (liberal)
2) 70%
3) exclude 3sd less Coverage
1 or 3, descriptive script

# ATLAS PATHS TO CHECK FOR RESLICING/COV CHECKS

# GET MATRIX OF VALID DATA FOR SESSIONS

# PATHS FOR EACH ART SESSION, OG RPTxt

1. Valid Scan MATRIX (face,reward1,reward2, rest)
2. rp
3. ART stuff
3. physio (when possible)
4. ATLAS checks

can conn support empty sessions (re valid scan matrix consistently lining up)
'''
