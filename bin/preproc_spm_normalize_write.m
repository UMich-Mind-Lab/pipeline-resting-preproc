function preproc_spm_normalize_write(varargin)

try
  % parse inputs
  p = inputParser;

  addParameter(p,'inDef','');
  addParameter(p,'inNii','');
  addParameter(p,'outFile','');
  addParameter(p,'spmConfig','config/spm_config.json');

  parse(p,varargin{:});

  %get number of volumes
  nVols = length(spm_vol(p.Results.inNii));

  %create matrix of functional volumes (e.g., 'func.nii,1'; 'func.nii,2'; ...)
  volumes = cell(nVols,1);
  for i = 1:nVols
    volumes{i,1} = strcat(p.Results.inNii,',',num2str(i));
  end

  %read spmConfig into workspace
  spmCfg = jsondecode(fileread(p.Results.spmConfig));
  spmCfg = spmCfg.(mfilename());

  mbatch{1}.spm.spatial.normalise.write.subj.def = {p.Results.inDef};
  mbatch{1}.spm.spatial.normalise.write.subj.resample = volumes;
  mbatch{1}.spm.spatial.normalise.write.woptions.bb = spmCfg.bb;
  mbatch{1}.spm.spatial.normalise.write.woptions.vox = spmCfg.vox;
  mbatch{1}.spm.spatial.normalise.write.woptions.interp = spmCfg.interp;
  mbatch{1}.spm.spatial.normalise.write.woptions.prefix = spmCfg.prefix;

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
