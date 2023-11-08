function experiment = loadExp()
% Load experiment object
cd /media/ben/'Extreme SSD'/analysisData/added_aligned/ % optional 
[file,path] = uigetfile('*.mat');
cd(path)
load(file)
end

