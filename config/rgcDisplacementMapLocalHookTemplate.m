function rgcDisplacementMapLocalHook
% OLFlickerSensitivityLocalHook - Configure things for working on OneLight projects.
%
% For use with the ToolboxToolbox.  If you copy this into your
% ToolboxToolbox localToolboxHooks directory (by default,
% ~/localToolboxHooks) and delete "LocalHooksTemplate" from the filename,
% this will get run when you execute tbUse({'OLFlickerSensitivityConfig'}) to set up for
% this project.  You then edit your local copy to match your local machine.
%
% The main thing that this does is define Matlab preferences that specify input and output
% directories.
%
% You will need to edit the project location and i/o directory locations
% to match what is true on your computer.

 
%% Define project
toolboxName = 'rgcDisplacementMap';
 
%% Clear out old preferences
if (ispref(toolboxName))
    rmpref(toolboxName);
end
 
%% Specify project location
toolboxBaseDir = tbLocateToolbox(toolboxName);

%% Set preferences for project output
setpref(toolboxName,'LocalDataPath',fullfile(toolboxBaseDir,'data')); % path to small file within the git repo 
end
