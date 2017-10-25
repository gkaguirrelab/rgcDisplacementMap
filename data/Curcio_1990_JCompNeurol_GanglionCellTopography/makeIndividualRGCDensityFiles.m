
tableFileName = 'DENSITY5_gc_resorted.csv';
outputFileNameStem = 'curcioRGCDensityPerSqMm_';

T=readtable(tableFileName);

SupportPosDeg=table2array(T(1,4:end));

uniqueSubjectNames = unique(table2array(T(2:end,1)));

meridianNames = {'NASAL','SUPERIOR','TEMPORAL','INFERIOR'};

for nn=1:length(uniqueSubjectNames)
    curcioRGCDensityPerSqMm = struct;
    curcioRGCDensityPerSqMm.supportPosDeg = SupportPosDeg;
    subjectName=uniqueSubjectNames{nn};
    for mm = 1:length(meridianNames)
        thisMeridianName = meridianNames{mm};
    rowIdx=find( strcmp(table2array(T(:,1)),subjectName) .* strcmp(table2array(T(:,2)),thisMeridianName) );
    curcioRGCDensityPerSqMm.(lower(thisMeridianName)) = table2array(T(rowIdx,4:end));
    end
    curcioRGCDensityPerSqMm.meta.subjectName = subjectName;
    curcioRGCDensityPerSqMm.meta.dataTableName = tableFileName;
    save([outputFileNameStem subjectName], 'curcioRGCDensityPerSqMm');
end
