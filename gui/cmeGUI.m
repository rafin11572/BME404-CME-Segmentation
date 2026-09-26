function app = cmeGUI()
%CMEGUI  Launch the BME 404 CME segmentation GUI.
%
%   >> cmeGUI
%
% Convenience wrapper so the app can be started with a single word during the
% demo.  All the interface code lives in CMEApp.m.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root,'src'));
addpath(here);
app = CMEApp();
if nargout == 0, clear app; end
end
