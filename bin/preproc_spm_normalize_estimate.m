function preproc_spm_normalize_estimate(varargin)

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

  mbatch{1}.spm.spatial.normalise.est.subj.vol = {p.Results.inStruct};
  mbatch{1}.spm.spatial.normalise.est.eoptions.biasreg = spmCfg.biasreg;
  mbatch{1}.spm.spatial.normalise.est.eoptions.biasfwhm = spmCfg.biasfwhm;
  mbatch{1}.spm.spatial.normalise.est.eoptions.tpm = {p.Results.TPM};
  mbatch{1}.spm.spatial.normalise.est.eoptions.affreg = spmCfg.affreg;
  mbatch{1}.spm.spatial.normalise.est.eoptions.reg = spmCfg.reg
  mbatch{1}.spm.spatial.normalise.est.eoptions.fwhm = spmCfg.fwhm;
  mbatch{1}.spm.spatial.normalise.est.eoptions.samp = spmCfg.samp;


  %save job to run in container
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
