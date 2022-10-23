function experiment = loadExp()
% Load experiment object
cd /media/ben/'Extreme SSD'/analysisData/ % optional 
[file,path] = uigetfile('*.mat');
cd(path)
load(file)
end

