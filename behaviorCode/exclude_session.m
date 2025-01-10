function TrialAnSt = exclude_session(TrialAnSt, exclude_crit)
% exclude_crit = {'MAW14' '06/01/19'; 'MAW14' '06/01/19'; 'MAW14' '06/01/19'};
for i_session = 1:size(exclude_crit,1)
    for i = 1:length(TrialAnSt)
        if ~isempty(TrialAnSt(i).(exclude_crit{i_session,1}))
            dates{i} = TrialAnSt(i).(exclude_crit{i_session,1})(1).mpc.StartDate;
        end
    end
    dateIdx = strcmp(exclude_crit{i_session,2},dates);
    TrialAnSt(dateIdx).(exclude_crit{i_session,1}) = [];
end
