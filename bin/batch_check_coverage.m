function batch_check_coverage(niiFile,maskStr,outDir)

try
  %p = inputParser;

  %addParameter(p,'niiFile','');
  %addParameter(p,'maskStr','');
  %addParameter(p,'outDir','');

  %parse(p,varargin);

  % split maskStr into cell array
  masks = split(maskStr,' ');

  %loop through masks
  for i = 1:length(masks)
    %get name of mask and append to outDir
    [~,maskName,~] = fileparts(masks{i});
    outFile = fullfile(outDir,strcat(maskName,'.txt'));
    if ~exist(outDir,'dir')
      mkdir(outDir);
    end
    %calculate coverage
    fprintf('nii: %s\nmask: %s\noutput: %s\n',niiFile,masks{i},outFile);
    check_coverage(niiFile,masks{i},outFile);
  end

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
