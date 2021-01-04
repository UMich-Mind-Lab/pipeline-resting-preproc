function preproc_spm_segment(varargin)


try
  % parse inputs
  p = inputParser;

  addParameter(p,'inStruct','');
  addParameter(p,'outFile','');
  addParameter(p,'TPM','');
  addParameter(p,'spmConfig','config/spm_config.json');
  parse(p,varargin{:});

  %read spmConfig into workspace
  spmCfg = jsondecode(fileread(p.Results.spmConfig));
  spmCfg = spmCfg.(mfilename());

  mbatch{1}.spm.spatial.preproc.channel.vols = {p.Results.inStruct};
  mbatch{1}.spm.spatial.preproc.channel.biasreg = 0;
  mbatch{1}.spm.spatial.preproc.channel.biasfwhm = Inf;
  mbatch{1}.spm.spatial.preproc.channel.write = [0 0];
  mbatch{1}.spm.spatial.preproc.tissue(1).tpm = {sprintf('%s,1',p.Results.TPM)};
  mbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
  mbatch{1}.spm.spatial.preproc.tissue(1).native = [0 0];
  mbatch{1}.spm.spatial.preproc.tissue(1).warped = [1 0];
  mbatch{1}.spm.spatial.preproc.tissue(2).tpm = {sprintf('%s,2',p.Results.TPM)};
  mbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
  mbatch{1}.spm.spatial.preproc.tissue(2).native = [0 0];
  mbatch{1}.spm.spatial.preproc.tissue(2).warped = [1 0];
  mbatch{1}.spm.spatial.preproc.tissue(3).tpm = {sprintf('%s,3',p.Results.TPM)};
  mbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
  mbatch{1}.spm.spatial.preproc.tissue(3).native = [0 0];
  mbatch{1}.spm.spatial.preproc.tissue(3).warped = [1 0];
  mbatch{1}.spm.spatial.preproc.tissue(4).tpm = {sprintf('%s,4',p.Results.TPM)};
  mbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
  mbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
  mbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
  mbatch{1}.spm.spatial.preproc.tissue(5).tpm = {sprintf('%s,5',p.Results.TPM)};
  mbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
  mbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
  mbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
  mbatch{1}.spm.spatial.preproc.tissue(6).tpm = {sprintf('%s,6',p.Results.TPM)};
  mbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
  mbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
  mbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
  mbatch{1}.spm.spatial.preproc.warp.mrf = spmCfg.mrf;
  mbatch{1}.spm.spatial.preproc.warp.cleanup = spmCfg.cleanup;
  mbatch{1}.spm.spatial.preproc.warp.reg = spmCfg.reg;
  mbatch{1}.spm.spatial.preproc.warp.affreg = spmCfg.affreg;
  mbatch{1}.spm.spatial.preproc.warp.fwhm = spmCfg.fwhm;
  mbatch{1}.spm.spatial.preproc.warp.samp = spmCfg.samp;
  mbatch{1}.spm.spatial.preproc.warp.write = spmCfg.write;
  mbatch{1}.spm.spatial.preproc.warp.vox = spmCfg.vox;
  mbatch{1}.spm.spatial.preproc.warp.bb = spmCfg.bb;

  %save job file
  save(p.Results.outFile,'mbatch');

catch ME
  fprintf('MATLAB code threw an exception:\n')
  fprintf('%s\n',ME.message);
  if length(ME.stack) ~= 0
    for i = 1:length(ME.stack)
      fprintf('File:%s\nName:%s\nLine:%d\n',ME.stack(i).file,...
        ME.stack(i).name,ME.stack(i).line);
    end
  end
end
