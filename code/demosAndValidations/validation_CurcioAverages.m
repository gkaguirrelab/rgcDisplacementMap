function validation_CurcioAverages(varargin)
% validation_CurcioAverages
%
% This routine compares the average cone and RGC densities reported in
% Christine Curcio's two 1990 J Comp Neurology papers to the averages that
% we compute from the raw individual data that she has provided.

%% Cone density functions

curcioReportedAverageConeDensity = [

calculatedAverageConeDensityFile = ...
fullfile([getpref('rgcDisplacementMap','LocalDataPath') , '/Curcio_1990_JCompNeurol_HumanPhotoreceptorTopography/curcioRawConeDensity_average.mat']), ...


end % function