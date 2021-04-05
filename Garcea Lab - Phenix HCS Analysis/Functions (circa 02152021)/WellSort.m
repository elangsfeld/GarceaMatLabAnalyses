function [Results_Sorted,WellNumber] = WellSort(Channels,Results)
for s = 1:size(Results,2)
    if Results(s).TotalNuclei > 0
        NoCellFilter(s,1) = 1;
    else
        NoCellFilter(s,1) = 0;
    end
end

Results_Filtered = Results(logical(NoCellFilter));

for i = 1:size(Results_Filtered,2)
    FieldsofView(i,1) = string(Results_Filtered(i).FieldofView);
    WellsMeasured(i,1) = string(Results_Filtered(i).Well);
end

WellsUnique = unique(WellsMeasured);
WellNumber = length(WellsUnique);

for w = 1:WellNumber
    Results_Sorted(w).Well = WellsUnique(w,1);
    NumFields(w,1) = sum(WellsMeasured == WellsUnique(w,1));
    Results_Sorted(w).NumFields = NumFields(w,1);
    for f = (sum(NumFields(1:(w-1)))+1):sum(NumFields(1:w,1))
        if f==(sum(NumFields(1:(w-1)))+1)
            Results_Sorted(w).TotalNumberofNuclei = Results_Filtered(f).TotalNuclei;
            Results_Sorted(w).MeanCh1(1:size(Results_Filtered(f).MeanIntensities,1),1) = Results_Filtered(f).MeanIntensities(:,1);
            Results_Sorted(w).MeanCh2(1:size(Results_Filtered(f).MeanIntensities,1),1) = Results_Filtered(f).MeanIntensities(:,2);
            if Channels >2, Results_Sorted(w).MeanCh3(1:size(Results_Filtered(f).MeanIntensities,1),1) = Results_Filtered(f).MeanIntensities(:,3); else end
            if Channels >3, Results_Sorted(w).MeanCh4(1:size(Results_Filtered(f).MeanIntensities,1),1) = Results_Filtered(f).MeanIntensities(:,4); else end
            Results_Sorted(w).SumCh1(1:size(Results_Filtered(f).SumIntensities,1),1) = Results_Filtered(f).SumIntensities(:,1);
            Results_Sorted(w).SumCh2(1:size(Results_Filtered(f).SumIntensities,1),1) = Results_Filtered(f).SumIntensities(:,2);
            if Channels >2, Results_Sorted(w).SumCh3(1:size(Results_Filtered(f).SumIntensities,1),1) = Results_Filtered(f).SumIntensities(:,3);else end
            if Channels >3, Results_Sorted(w).SumCh4(1:size(Results_Filtered(f).SumIntensities,1),1) = Results_Filtered(f).SumIntensities(:,4); else end            
        else
            Results_Sorted(w).TotalNumberofNuclei = Results_Filtered(f).TotalNuclei + Results_Sorted(w).TotalNumberofNuclei;
            Results_Sorted(w).TotalNuclei(Results_Filtered(f).FieldofView,1) = Results_Filtered(f).TotalNuclei;
            Results_Sorted(w).MeanCh1(end+1:end+size(Results_Filtered(f).MeanIntensities,1),1) = Results_Filtered(f).MeanIntensities(:,1);
            Results_Sorted(w).MeanCh2(end+1:end+size(Results_Filtered(f).MeanIntensities,1),1) = Results_Filtered(f).MeanIntensities(:,2);
            if Channels >2, Results_Sorted(w).MeanCh3(end+1:end+size(Results_Filtered(f).MeanIntensities,1),1) = Results_Filtered(f).MeanIntensities(:,3); else end
            if Channels >3, Results_Sorted(w).MeanCh4(end+1:end+size(Results_Filtered(f).MeanIntensities,1),1) = Results_Filtered(f).MeanIntensities(:,4); else end
            Results_Sorted(w).SumCh1(end+1:end+size(Results_Filtered(f).SumIntensities,1),1) = Results_Filtered(f).SumIntensities(:,1);
            Results_Sorted(w).SumCh2(end+1:end+size(Results_Filtered(f).SumIntensities,1),1) = Results_Filtered(f).SumIntensities(:,2);
            if Channels >2, Results_Sorted(w).SumCh3(end+1:end+size(Results_Filtered(f).SumIntensities,1),1) = Results_Filtered(f).SumIntensities(:,3); else end
            if Channels >3, Results_Sorted(w).SumCh4(end+1:end+size(Results_Filtered(f).SumIntensities,1),1) = Results_Filtered(f).SumIntensities(:,4); else end
        end
    end
end
end