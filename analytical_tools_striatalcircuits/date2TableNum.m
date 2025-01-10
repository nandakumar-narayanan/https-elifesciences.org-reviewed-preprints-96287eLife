
%returns an index for the day??

function dateNum = date2TableNum(animal, date, TrialAnSt)
if isfield(TrialAnSt,animal)
    animals = fieldnames(TrialAnSt);
    for animalIdx = 1:length(animals)
        animal_id = char(animals(animalIdx));
        for dayIdx = 1:size(TrialAnSt,2)
            if ~isempty(TrialAnSt(dayIdx).(animal))
                if strcmp(date,TrialAnSt(dayIdx).(animal)(1).mpc.StartDate)
                    dateNum = dayIdx;
                end
            end
        end
    end
else
    return
end
end