function loadMat()
[file,path] = uigetfile('*.mat');
cd(path)
load(file)
end
