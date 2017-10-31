
tableFileNames = {'DENSITY8_cones_resorted_computedAverage.xlsx','curcioReportedConeAverage.xlsx'};
outputFileNameStem = 'curcioRawConeDensity_';

for tt = 1:length(tableFileNames)
    T=readtable(tableFileNames{tt});
    
    support=table2array(T(1,4:end));
    
    uniqueSubjectNames = unique(table2array(T(2:end,1)));
    
    meridianNames = {'NASAL','SUPERIOR','TEMPORAL','INFERIOR'};
    
    for nn=1:length(uniqueSubjectNames)
        curcioRawConeDensity = struct;
        curcioRawConeDensity.support = support;
        subjectName=uniqueSubjectNames{nn};
        for mm = 1:length(meridianNames)
            thisMeridianName = meridianNames{mm};
            rowIdx=find( strcmp(table2array(T(:,1)),subjectName) .* strcmp(table2array(T(:,2)),thisMeridianName) );
            curcioRawConeDensity.(lower(thisMeridianName)) = table2array(T(rowIdx,4:end));
        end
        curcioRawConeDensity.meta.subjectName = upper(subjectName);
        curcioRawConeDensity.meta.dataTableName = tableFileNames{tt};
        curcioRawConeDensity.meta.densityUnits = T.Units{rowIdx};
        curcioRawConeDensity.meta.supportUnits = T.Units{1};
        save([outputFileNameStem upper(subjectName)], 'curcioRawConeDensity');
    end % loop over subjects
    
end % loop over tables