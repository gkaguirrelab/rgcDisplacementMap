
tableFileName = 'DENSITY8_cones_resorted.csv';
outputFileNameStem = 'curcioConeDensityPerSqMm_';

T=readtable(tableFileName);

SupportPosDeg=table2array(T(1,4:end));

uniqueSubjectNames = unique(table2array(T(2:end,1)));

meridianNames = {'NASAL','SUPERIOR','TEMPORAL','INFERIOR'};

for nn=1:length(uniqueSubjectNames)
    curcioConeDensityPerSqMm = struct;
    curcioConeDensityPerSqMm.supportPosDeg = SupportPosDeg;
    subjectName=uniqueSubjectNames{nn};
    for mm = 1:length(meridianNames)
        thisMeridianName = meridianNames{mm};
    rowIdx=find( strcmp(table2array(T(:,1)),subjectName) .* strcmp(table2array(T(:,2)),thisMeridianName) );
    curcioConeDensityPerSqMm.(lower(thisMeridianName)) = table2array(T(rowIdx,4:end));
    end
    curcioConeDensityPerSqMm.meta.subjectName = subjectName;
    curcioConeDensityPerSqMm.meta.dataTableName = tableFileName;
    save([outputFileNameStem subjectName], 'curcioConeDensityPerSqMm');
end
