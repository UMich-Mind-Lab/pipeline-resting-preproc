function preproc_spm_realignunwarp(varargin)

try
  % parse inputs
  p = inputParser;

  validChar = @(x) ischar(x) && (x>0);
  addParameter(p,'inFunc',validChar);
  addParameter(p,'inVDM',validChar)
  addParameter(p,'outFile',validChar);
  parse(p,varargin{:});

  %get n Volumes for functional image
  nVols = length(spm_vol(p.Results.inFunc));

  %create matrix of functional volumes (e.g., 'func.nii,1'; 'func.nii,2'; ...)
  volumes = cell(nVols,1);
  for i = 1:nVols
    volumes{i,1} = strcat(p.Results.inFunc,',',num2str(i));
  end

  mbatch{1}.spm.spatial.realignunwarp.data.scans = volumes;
  mbatch{1}.spm.spatial.realignunwarp.data.pmscan = {strcat(p.Results.inVDM,',1')};
  mbatch{1}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
  mbatch{1}.spm.spatial.realignunwarp.eoptions.sep = 4;
  mbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
  mbatch{1}.spm.spatial.realignunwarp.eoptions.rtm = 0;
  mbatch{1}.spm.spatial.realignunwarp.eoptions.einterp = 2;
  mbatch{1}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
  mbatch{1}.spm.spatial.realignunwarp.eoptions.weight = '';
  mbatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
  mbatch{1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
  mbatch{1}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
  mbatch{1}.spm.spatial.realignunwarp.uweoptions.jm = 0;
  mbatch{1}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
  mbatch{1}.spm.spatial.realignunwarp.uweoptions.sot = [];
  mbatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
  mbatch{1}.spm.spatial.realignunwarp.uweoptions.rem = 1;
  mbatch{1}.spm.spatial.realignunwarp.uweoptions.noi = 5;
  mbatch{1}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
  mbatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
  mbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
  mbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
  mbatch{1}.spm.spatial.realignunwarp.uwroptions.mask = 1;
  mbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';

  %prep job to run through container
  save(p.Results.outFile,'mbatch')

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
