
tableFileName = 'DENSITY8_cones_resorted_computedAverage.xlsx';
outputFileNameStem = 'curcioRawConeDensity_';

T=readtable(tableFileName);

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
    curcioRawConeDensity.meta.subjectName = subjectName;
    curcioRawConeDensity.meta.dataTableName = tableFileName;
            curcioRawConeDensity.meta.densityUnits = T.Units{rowIdx};
        curcioRawConeDensity.meta.supportUnits = T.Units{1};
    save([outputFileNameStem subjectName], 'curcioRawConeDensity');
end
