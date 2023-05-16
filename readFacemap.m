function [keypoint_out] = readFacemap(filename,keypoint)
% This function takes the data from facemap h5 file of the keypoint
% positions and returns a structure with the x,y coordinates
keypoint_out.x = h5read(filename,['/Facemap/',keypoint,'/x']);
keypoint_out.y = h5read(filename,['/Facemap/',keypoint,'/y']);
end