function L1_spm_modelcon(varargin)

try
  % parse inputs
  p = inputParser;

  validChar = @(x) ischar(x) && (x>0);
  addParameter(p,'SPM',validChar);
  addParameter(p,'task',validChar);
  addParameter(p,'outFile',validChar)

  parse(p,varargin{:});

  %load SPM file into workspace
  load(p.Results.SPM);

  %extract condition names
  allConditions = {};
  nSess = length(SPM.Sess);
  for i = 1:nSess
    allConditions(i) = {vertcat(SPM.Sess(i).U.name)};
  end

  if strcmp(p.Results.task,'reward')
    %% define Anticipation contrasts
    config.con{1}.name = 'WinAnt>LossAnt';
    config.con{1}.conditions = {'WinBigAnticipation','WinSmallAnticipation','LoseBigAnticipation','LoseSmallAnticipation';
                                'WinBigAnticipation','WinSmallAnticipation','LoseBigAnticipation','LoseSmallAnticipation'};
    config.con{1}.weights = [.25,.25,-.25,-.25;
                             .25,.25,-.25,-.25];
    config.con{2}.name = 'WinAnt>NeutAnt';
    config.con{2}.conditions = {'WinBigAnticipation','WinSmallAnticipation', 'NeutralAnticipation';
                                'WinBigAnticipation','WinSmallAnticipation', 'NeutralAnticipation'};
    config.con{2}.weights = [.25,.25,-.5;
                             .25,.25,-.5];
    config.con{3}.name = 'LossAnt>NeutAnt';
	config.con{3}.conditions = {'LoseBigAnticipation','LoseSmallAnticipation', 'NeutralAnticipation';
	                            'LoseBigAnticipation','LoseSmallAnticipation', 'NeutralAnticipation'};	
    config.con{3}.weights = [.25,.25,-.5;
	                         .25,.25,-.5];
	config.con{4}.name = 'BigWinAnt>NeutAnt';
	config.con{4}.conditions = {'WinBigAnticipation','NeutralAnticipation';
	                            'WinBigAnticipation','NeutralAnticipation'};
    config.con{4}.weights = [.5,-.5;
                             .5,-.5];
    config.con{5}.name = 'SmallWinAnt>NeutAnt';
    config.con{5}.conditions = {'WinSmallAnticipation','NeutralAnticipation';
                                'WinSmallAnticipation','NeutralAnticipation'};
    config.con{5}.weights = [.5,-.5;
                             .5,-.5];
    config.con{6}.name = 'BigLossAnt>NeutAnt'
    config.con{6}.conditions = {'LoseBigAnticipation','NeutralAnticipation';
                                'LoseBigAnticipation','NeutralAnticipation'};
    config.con{6}.weights = [.5,-.5;
                             .5,-.5];
    config.con{7}.name = 'SmallLossAnt>NeutAnt'
    config.con{7}.conditions = {'LoseSmallAnticipation','NeutralAnticipation';
                                'LoseSmallAnticipation','NeutralAnticipation'};
    config.con{7}.weights = [.5,-.5;
                             .5,-.5];
	config.con{8}.name = 'BigWinAnt>SmallWinAnt'
	config.con{8}.conditions = {'WinBigAnticipation','WinSmallAnticipation';
                                'WinBigAnticipation','WinSmallAnticipation'};
    config.con{8}.weights = [.5,-.5;
                             .5,-.5];
	config.con{9}.name = 'BigLossAnt>SmallLossAnt'
    config.con{9}.conditions = {'LoseBigAnticipation','LoseSmallAnticipation';
                                'LoseBigAnticipation','LoseSmallAnticipation'};
    config.con{9}.weights = [.5,-.5;
                             .5,-.5];
	%%contrasts that compared run 1 vs run 2
    % config.con{2}.name = 'WinFeedbackPos1>WinFeedbackPos2';
    % config.con{2}.conditions = {'WinBigFeedbackPositive','WinSmallFeedbackPositive';
                                % 'WinBigFeedbackPositive','WinSmallFeedbackPositive'};
    % config.con{2}.weights = [.5,.5;
                             % -.5,-.5];

  end

  %% generate matlabbatch
  % loop through each contrast in config
  for i = 1:length(config.con)
      % verify that all conditions specified exist
      hasAllConditions=true;
      for j = 1:nSess
        hasAllConditions = hasAllConditions & all(ismember(config.con{i}.conditions(j,:),allConditions{j}));
      end
      if hasAllConditions
        % initialize weights array
        w = zeros(1,length(vertcat(allConditions{:})));
        % loop through sessions
        for j = 1:nSess
          %extract conditions for current session
          c=config.con{i}.conditions(j,:);
          %if not first session, then all condition indices need to be shifted
          %based on how many conditions existed in the prior session
          nShift=0;
          if j > 1
            nShift=length(vertcat(allConditions{1:j-1}));
          end
          %now loop through conditions for current session
          for k = 1:length(c)
            %get index for condition within session condition list
            idx = find(contains(allConditions{j},c{k})) + nShift;
            %add appropriate value to weights array
            w(idx) = config.con{i}.weights(j,k);
          end
        end
        %populate mbatch for current contrast
        mbatch{1}.spm.stats.con.consess{i}.tcon.name = config.con{i}.name;
        mbatch{1}.spm.stats.con.consess{i}.tcon.weights = w;
        mbatch{1}.spm.stats.con.consess{i}.tcon.sessrep = 'none';
      else
          fprintf('uh oh spaghetti-o! No %s contrast for you!\n',config.con{i}.name);
      end
  end

  %populate rest of batch info
  mbatch{1}.spm.stats.con.spmmat = {p.Results.SPM};
  mbatch{1}.spm.stats.con.delete = 1;

  %% save matlabbatch to run in container
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
