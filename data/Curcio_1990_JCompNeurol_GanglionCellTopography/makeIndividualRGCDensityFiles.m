
tableFileName = 'DENSITY5_gc_resorted_computedAverage.xls';
outputFileNameStem = 'curcioRawRGCDensity_';

T=readtable(tableFileName);

SupportPosDeg=table2array(T(1,4:end));

uniqueSubjectNames = unique(table2array(T(2:end,1)));

meridianNames = {'NASAL','SUPERIOR','TEMPORAL','INFERIOR'};

for nn=1:length(uniqueSubjectNames)
    curcioRawRGCDensity = struct;
    curcioRawRGCDensity.support = SupportPosDeg;
    subjectName=uniqueSubjectNames{nn};
    for mm = 1:length(meridianNames)
        thisMeridianName = meridianNames{mm};
    rowIdx=find( strcmp(table2array(T(:,1)),subjectName) .* strcmp(table2array(T(:,2)),thisMeridianName) );
    curcioRawRGCDensity.(lower(thisMeridianName)) = table2array(T(rowIdx,4:end));
    end
    curcioRawRGCDensity.meta.subjectName = subjectName;
    curcioRawRGCDensity.meta.dataTableName = tableFileName;
    curcioRawRGCDensity.meta.densityUnits = 'counts/mm2';
    curcioRawRGCDensity.meta.supportUnits = 'deg';
    save([outputFileNameStem subjectName], 'curcioRawRGCDensity');
end
