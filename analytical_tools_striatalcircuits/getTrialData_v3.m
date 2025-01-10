%% Code built to pull switch data from MedPC files
%% Adjusted from Youngcho's getTrialData and getTrialData_Bisection
%% Austin Bruce, Tomas Lence, and Matt Weber
%% Clarifications added by Kumar Narayanan 9/2024

%% Changes to make - notes go here

%% Med Associates Code
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
\   R	    Opto laser ; check line 84-87
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
function TrialAnSt = getTrialData_v3(mpcParsed)

uniqueSbj = unique({mpcParsed.Subject});
nmpcType = min(cellfun(@length, {mpcParsed.H}), cellfun(@length, {mpcParsed.S}))> 0; % check if any trials

%it's CR5 with the error ###
for i_animal = 1:size(uniqueSbj,2)
    sbjIdx = strcmp(uniqueSbj(i_animal),{mpcParsed.Subject}) &  nmpcType;
    sbjLineIdx = find(sbjIdx);
    numOfDays = sum(sbjIdx);
    animals = char(uniqueSbj(i_animal));
    
    for i_day = 1:numOfDays
        lineIdx =sbjLineIdx(i_day);
        if contains(mpcParsed(lineIdx).MSN,'18L6R')
            type = 1;
        elseif contains(mpcParsed(lineIdx).MSN,'6L18R')
            type = 0;
        else
            type = NaN;
        end
        
        timeITI = mpcParsed(lineIdx).U';
        trialStart = mpcParsed(lineIdx).H'; % 1
        trialEnd = mpcParsed(lineIdx).S'; % 2
        trialStart=trialStart(1:length(trialEnd));
        trialType = mpcParsed(lineIdx).T'; % the programmed duration for the trial
        trialType=trialType(1:length(trialEnd));
        reward = mpcParsed(lineIdx).W';
%         Note counter-intuitive coding; all animals verified by direct
%         observation prior to session; example laserTest.MOV in /MPChalo
%         This note added after questions from Paton Lab 8/2024
          if ismember(uniqueSbj(i_animal), {'CR2', 'CR3', 'CL3', 'SL3', 'SL6', 'SL8', 'SL9', 'SR1', 'SR2' })
                opto = mpcParsed(lineIdx).R'; opto = opto==2; % with the new arduino protocol to use after the lab moved, this changed.
            else
                opto = mpcParsed(lineIdx).R'; opto = opto==1; % 1 means laser is off. 2 is on/ 
            end
            
            
            rpInLeft = mpcParsed(lineIdx).P'; % responding on left NP
            rpOutLeft = mpcParsed(lineIdx).Q';
            rpInRight = mpcParsed(lineIdx).N'; % responding on right NP
            rpOutRight = mpcParsed(lineIdx).O';
            rpInBack = mpcParsed(lineIdx).F'; % responding on right NP
            %         rpOutBack = mpcParsed(lineIdx).G'; % correct: being used for NP out and side choice...
            %         fzIn = mpcParsed(lineIdx).Z';
            %         fzOut = mpcParsed(lineIdx).Y';
            
            trialOutcomeInfo = mpcParsed(lineIdx).L'; % 4 = correct trial, 3 = incorrect trial, 2 = correct correction trial, 1 = incorrect correction trial
            if type == 0    % 0 indicates that a short latency trial is rewarded at the left nose poke
                rsp_in_short = rpInLeft;
                rsp_out_short = rpOutLeft;
            elseif type == 1
                rsp_in_short = rpInRight;
                rsp_out_short =rpOutRight;
            end
            
%             if numel(rsp_in_short) == numel(rsp_out_short)
%                 rsp_short_duration = rsp_out_short - rsp_in_short;
%             elseif numel(rsp_in_short) == numel(rsp_out_short)+1 && rsp_in_short(end) > rsp_out_short(end)
%                 rsp_short_duration = rsp_out_short - rsp_in_short(1:end-1);
%             end
%             
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
                
                %% response time within trial -- check here release time is incorrect
                trial(j).leftRspTimeTrial = rpInLeft(rpInLeft>curTS & rpInLeft<=curTE) - curTS;
                trial(j).leftRelTimeTrial = rpOutLeft(rpOutLeft>curTS & rpOutLeft<=curTE+10) - curTS; % one second used as a buffer time
                % select those under length of responses
                if(~isempty(trial(j).leftRelTimeTrial) & (length(trial(j).leftRelTimeTrial)>length(trial(j).leftRspTimeTrial)))
                    trial(j).leftRelTimeTrial = trial(j).leftRelTimeTrial(1:length(trial(j).leftRspTimeTrial));
                end
                trial(j).rightRspTimeTrial = rpInRight(rpInRight>curTS & rpInRight<=curTE) - curTS;
                trial(j).rightRelTimeTrial = rpOutRight(rpOutRight>curTS & rpOutRight<=curTE+10) - curTS; % one second used as a buffer time
                % select those under length of responses
                if (~isempty(trial(j).rightRelTimeTrial) & (length(trial(j).rightRelTimeTrial)>length(trial(j).rightRspTimeTrial)))
                    trial(j).rightRelTimeTrial = trial(j).rightRelTimeTrial(1:length(trial(j).rightRspTimeTrial));
                end
                
                if type == 0    % 0 indicates that a short latency trial is rewarded at the left nose poke
                    short_rsp = trial(j).leftRspTimeTrial;
                    short_rel = trial(j).leftRelTimeTrial;
                    long_rsp = trial(j).rightRspTimeTrial;
                    long_rel = trial(j).rightRelTimeTrial;
                elseif type == 1
                    short_rsp = trial(j).rightRspTimeTrial;
                    short_rel = trial(j).rightRelTimeTrial;
                    long_rsp = trial(j).leftRspTimeTrial;
                    long_rel = trial(j).leftRelTimeTrial;
                else
                    error('')
                end
                switch_arrival = [];
                switch_depart = [];
                if ~isempty(short_rel) && ~isempty(long_rsp) && any(long_rsp > min(short_rel)) && trialType(j) == 18000
                    switch_arrival = min(long_rsp(long_rsp > min(short_rel)));
                    switch_depart = max(short_rel(short_rel < switch_arrival));
                end
                % Include the portion below if you want to include short trial switches
                %             if ~isempty(short_rel) && ~isempty(long_rsp) && any(long_rsp > min(short_rel)) && trialType(j) == 6000
                %                 switch_arrival = min(long_rsp(long_rsp > min(short_rel)));
                %                 switch_depart = max(short_rel(short_rel < switch_arrival));
                %             end
                
                trial(j).SwitchArrival = switch_arrival;
                trial(j).SwitchDepart = switch_depart;
                trial(j).ShortRsp = short_rsp;
                trial(j).ShortRel = short_rel;
                trial(j).LongRsp = long_rsp;
                trial(j).LongRel = long_rel;
                
                firstPoke = [min(trial(j).rightRspTimeTrial) min(trial(j).leftRspTimeTrial)];
                trial(j).firstPoke = min(firstPoke);
%                 %% time from ITI end to trial initiation poke
%                 x = rpInBack - trialStart(j);
%                 trial(j).initiationRT = min(x(x>=0)); % timestamp of first response after the trial was initiated
                
                %% reaction time: time from signal offset to first response
                startPoke = min(rpInBack(rpInBack >= trialStart(j))); % timestamp of first back poke after trial start
                durationElapsed = startPoke+(trialType(j)./1000); % time that the duration elapsed
                potentialRTs = [min(rpInLeft(rpInLeft>=durationElapsed)) min(rpInRight(rpInRight>=durationElapsed))]; % finds first responses after duration elapsed on left and right nosepokes
                trial(j).RT = min(potentialRTs) - durationElapsed;
                
                %% reward entries/exits relative to trial end
                %             trial(j).rewardEntryTimes = fzIn - curTE;
                %             trial(j).rewardExitTimes = fzIn - curTE;
                trial(j).reward = reward(reward>curTS & reward<=curTE+0.2) - curTS;
            end
            %trial(1).rsp_short_duration = rsp_short_duration;
            trial(1).mpc = mpcParsed(lineIdx); % just saves everything
            TrialAnSt(i_day).(animals) = trial; %
            %TrialAnSt(i_day).(animals)([TrialAnSt(i_day).(animals).trialStart] > 3.6e+03) = [];
        end
        
    end
    
end