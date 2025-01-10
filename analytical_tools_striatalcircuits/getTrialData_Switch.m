%% Code built to pull switch data from MedPC files
%% Austin Bruce, Tomas Lence, and Matt Weber

%% Changes to make - notes go here

%% Not sure what this is for - remnant of Ben's 
%change choice configuration list to 12 and 21 for left solid / right blink and left blink / right solid
%Add G as a variable and SET G(C(^trialCounter)) = 1 for left choices and SET G(C(^trialCounter)) = 2; for right pokes make sure everything is sealed
% S8, \ correction trial
%     0.01": 
%         IF (C(^goodTrial) >= 1)  [@LOGDATA;@IGNOREDATA] \ logging correct/incorrect/correction like this will make it easier to integrate with matlab (i.e., goodtrialsindex=L>0)
% 	   @LOGDATA:    SET L(C(^trialCounter)) = 2; ---> SX; \ log that this non-correction trial was correct
% 	   @IGNOREDATA: SET L(C(^trialCounter)) = 0; ---> SX; !!!\ log that this was a correction trial. 
% 
%         ADD C(^trialCounter); 
%         SET C(^CorrectionTrial) = 1 ---> S9 

%% following are data format from MEDPC output file
%{
From Bisection task

\   A       Probe trial selection array
\   B       Reward counter
\   D
\   E       Interval Record Array
\   F
\   G       center nosepoke in array
\   H       Trial Start time record array
\   I       Variable Interval arrary
\   K       List of random ITI range
\   L       Log good trial accuracy data
\   M       Correction Trial tracker

\   N	    RIGHT NP response time record array
\   O       RIGHT NP release time record array

\   P       LEFT NP Response Time (Beam Break) record array
\   Q       LEFT NP Response Time (Beam Release) record array
\   R	    Opto laser on (1) or off (2)
\   S       Trial End time record array
\   T       Probe trial record array
\   U       ITI record array
\   V	    Probe length
\   W       Pellet dispense time record array
\   X       milli second Timer for recoring events
\   Y       Reward zone out record array
\   Z       Reward zone in record array
\-----------------------------------------
\----------------------------------------------------------
%}

%% Starts to pull values from MedPC files
function TrialAnSt = getTrialData(mpcParsed, typecell)
uniqueSbj = unique({mpcParsed.Subject});
ompcType = cellfun(@length, {mpcParsed.Q}) == 1 & cellfun(@length, {mpcParsed.Y}) == 0 & cellfun(@length, {mpcParsed.C}) == 0 & cellfun(@length, {mpcParsed.S}) > 0; % old MPC with trials
nmpcType = ~ompcType & min(cellfun(@length, {mpcParsed.H}), cellfun(@length, {mpcParsed.S}))> 0; % check if any trials
for i_animal = 1:size(uniqueSbj,2)
    sbjIdx = strcmp(uniqueSbj(i_animal),{mpcParsed.Subject}) & (ompcType | nmpcType);
    sbjLineIdx = find(sbjIdx);
    numOfDays = sum(sbjIdx);
    animals = char(uniqueSbj(i_animal));  
    animalsI = strcmp(animals,typecell(:,1));
    animalsIdx = find(animalsI);
    if ~isempty(animalsIdx)
    type = typecell{strcmp(animals,typecell(:,1)),2};
    for i_day = 1:numOfDays
        lineIdx =sbjLineIdx(i_day);

        timeITI = mpcParsed(lineIdx).U';
        trialStart = mpcParsed(lineIdx).H'; % 1
        trialEnd = mpcParsed(lineIdx).S'; % 2
        trialStart=trialStart(1:length(trialEnd));
        trialType = mpcParsed(lineIdx).T'; % the programmed duration for the trial
        trialType=trialType(1:length(trialEnd));
        reward = mpcParsed(lineIdx).W';
        opto = mpcParsed(lineIdx).R'; 
        opto = opto==1; % 1 means laser is on. 2 is off
            
        rpInLeft = mpcParsed(lineIdx).P'; % responding on left NP
        rpOutLeft = mpcParsed(lineIdx).Q';
        rpInRight = mpcParsed(lineIdx).N'; % responding on right NP
        rpOutRight = mpcParsed(lineIdx).O';
        rpInBack = mpcParsed(lineIdx).F'; % responding on right NP
        rpOutBack = mpcParsed(lineIdx).G'; % correct: being used for NP out and side choice...
        fzIn = mpcParsed(lineIdx).Z';
        fzOut = mpcParsed(lineIdx).Y';

        trialOutcomeInfo = mpcParsed(lineIdx).L'; % 4 = correct trial, 3 = incorrect trial, 2 = correct correction trial, 1 = incorrect correction trial
            
        trial = struct;
        trialNum = min(length(trialStart), length(trialEnd));
        for j = 1:trialNum
            %% trial start/end time, duration, and ITI
            curTS = trialStart(j); % H 
            curTE = trialEnd(j);
            trial(j).trialStart = curTS;
            trial(j).trialEnd = curTE;
            trial(j).trialDuration = curTE - curTS;
            trial(j).ITI = timeITI(j);
            
            %% real trial start time (when mouse iniated)
            x = rpInBack-curTS;
            trial(j).initiationRT = min(x(x>=0)); % timestamp of first response after the trial was initiated
            curTS =  trial(j).initiationRT +  curTS; % add the reaction time to trial start (when back light turned on), get the real trial start. 
            trial(j).realTrialStart = curTS;
            
            %% programmed duration for the trial
            trial(j).programmedDuration = trialType(j);
            
            %% opto
             if numel(opto)~=0; trial(j).opto = opto(j); end
             
            %% outcome info
            trial(j).outcome = trialOutcomeInfo(j);
            
            %% response time within trial
            trial(j).leftRspTimeTrial = rpInLeft(rpInLeft>curTS & rpInLeft<=curTE) - curTS;
            trial(j).leftRelTimeTrial = rpOutLeft(rpOutLeft>curTS & rpOutLeft<=curTE) - curTS;
            trial(j).rightRspTimeTrial = rpInRight(rpInRight>curTS & rpInRight<=curTE) - curTS;
            trial(j).rightRelTimeTrial = rpOutRight(rpOutRight>curTS & rpOutRight<=curTE) - curTS;
            
            %% Builds the Switch Departure and Switch Arrival Times based on behavioral training paradigm
            % 1. Start at short and then switch / 2. Start at short and then switch, back to short / 3. Multiple switches
            % 4. Start at long, move to short and then switch / 5. Start and stay at long / 6. Start and stay at short  
            if type == 0    % 0 indicates that a short latency trial is rewarded at the left nose poke 
                if ~isempty(trial(j).leftRelTimeTrial) && ~isempty(trial(j).rightRspTimeTrial) 
                    trial(j).SwitchArrival = min(trial(j).rightRspTimeTrial(trial(j).rightRspTimeTrial > min(trial(j).leftRelTimeTrial)));
                        if ~isempty(trial(j).SwitchArrival) 
                            trial(j).SwitchDepart = max(trial(j).leftRelTimeTrial(trial(j).leftRelTimeTrial < trial(j).SwitchArrival));
                        else trial(j).SwitchArrival=[]; 
                        end
                elseif isempty(trial(j).leftRelTimeTrial) && isempty(trial(j).rightRspTimeTrial)
                    trial(j).SwitchArrival = [];
                    trial(j).SwitchDepart = [];
                elseif ~isempty(trial(j).leftRelTimeTrial) && isempty(trial(j).rightRspTimeTrial)
                    trial(j).SwitchArrival = [];
                    trial(j).SwitchDepart = [];
                % Use this for the standard switch strategy - silence if using a different strategy
                elseif isempty(trial(j).leftRelTimeTrial) && ~isempty(trial(j).rightRspTimeTrial)
                    trial(j).SwitchArrival = [];
                    trial(j).SwitchDepart = [];
                % Use this for the no short response strategy - cautiously interpret this data
%                 elseif isempty(trial(j).leftRelTimeTrial) && ~isempty(trial(j).rightRspTimeTrial)
%                     trial(j).SwitchArrival = [];
%                     trial(j).SwitchDepart = min(trial(j).rightRspTimeTrial(trial(j).rightRspTimeTrial < trial(j).programmedDuration/1000));
                elseif isempty(trial(j).leftRelTimeTrial) && any(trial(j).rightRspTimeTrial > trial(j).programmedDuration/1000)
                    trial(j).SwitchArrival = [];
                    trial(j).SwitchDepart = [];
                else
                    trial(j).SwitchArrival = [];
                    trial(j).SwitchDepart = [];
                end
                trial(j).ShortRsp = trial(j).leftRspTimeTrial;
                trial(j).LongRsp = trial(j).rightRspTimeTrial;
                
            % 1. Start at short and then switch / 2. Start at short and then switch, back to short / 3. Multiple switches
            % 4. Start at long, move to short and then switch / 5. Start and stay at long / 6. Start and stay at short  
            elseif type == 1    % 1 indicates that a short latency trial is rewarded at the right nose poke             
                if ~isempty(trial(j).rightRelTimeTrial) && ~isempty(trial(j).leftRspTimeTrial) 
                    trial(j).SwitchArrival = min(trial(j).leftRspTimeTrial(trial(j).leftRspTimeTrial > min(trial(j).rightRelTimeTrial)));
                        if ~isempty(trial(j).SwitchArrival) 
                            trial(j).SwitchDepart = max(trial(j).rightRelTimeTrial(trial(j).rightRelTimeTrial < trial(j).SwitchArrival));
                        else trial(j).SwitchArrival=[]; 
                        end
                elseif isempty(trial(j).rightRelTimeTrial) && isempty(trial(j).leftRspTimeTrial)
                    trial(j).SwitchArrival = [];
                    trial(j).SwitchDepart = [];
                elseif ~isempty(trial(j).rightRelTimeTrial) && isempty(trial(j).leftRspTimeTrial)
                    trial(j).SwitchArrival = [];
                    trial(j).SwitchDepart = [];
                % Use this for the standard switch strategy - silence if using a different strategy
                elseif isempty(trial(j).rightRelTimeTrial) && ~isempty(trial(j).leftRspTimeTrial)
                    trial(j).SwitchArrival = [];
                    trial(j).SwitchDepart = [];
                % Use this for the no short response strategy - cautiously interpret this data
%                 elseif isempty(trial(j).rightRelTimeTrial) && ~isempty(trial(j).leftRspTimeTrial)
%                     trial(j).SwitchArrival = [];
%                     trial(j).SwitchDepart = min(trial(j).leftRspTimeTrial(trial(j).leftRspTimeTrial < trial(j).programmedDuration/1000);
                elseif isempty(trial(j).rightRelTimeTrial) && any(trial(j).leftRspTimeTrial > trial(j).programmedDuration/1000)
                    trial(j).SwitchArrival = [];
                    trial(j).SwitchDepart = [];
                else
                    trial(j).SwitchArrival = [];
                    trial(j).SwitchDepart = [];
                end
                trial(j).ShortRsp = trial(j).rightRspTimeTrial;
                trial(j).LongRsp = trial(j).leftRspTimeTrial;
            else
                error('')
            end               
                
            %% response times relative to trial start - **might be able to take this out**
            trial(j).leftResponseTimes = rpInLeft - curTS;
            trial(j).rightResponseTimes = rpInRight - curTS;
            trial(j).backResponseTimes = rpInBack - curTS;
            % release times relative to trial start
            trial(j).leftReleaseTimes = rpOutLeft - curTS;
            trial(j).rightReleaseTimes = rpOutRight - curTS;
            trial(j).backReleaseTimes = rpOutBack - curTS;
            
            %% time from ITI end to trial initiation poke
            x = rpInBack - trialStart(j);
            trial(j).initiationRT = min(x(x>=0)); % timestamp of first response after the trial was initiated
            
            %% reaction time: time from signal offset to first response
            startPoke = min(rpInBack(rpInBack >= trialStart(j))); % timestamp of first back poke after trial start
            durationElapsed = startPoke+(trialType(j)./1000); % time that the duration elapsed
            potentialRTs = [min(rpInLeft(rpInLeft>=durationElapsed)) min(rpInRight(rpInRight>=durationElapsed))]; % finds first responses after duration elapsed on left and right nosepokes
            trial(j).RT = min(potentialRTs) - durationElapsed; 

            %% reward entries/exits relative to trial end
            trial(j).rewardEntryTimes = fzIn - curTE;
            trial(j).rewardExitTimes = fzIn - curTE;
            trial(j).reward = reward;
            
        end
        trial(1).mpc = mpcParsed(lineIdx); % just saves everything
        TrialAnSt(i_day).(animals) = trial; %
    end
    else
    end

end