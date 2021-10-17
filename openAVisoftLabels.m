function [data] = openAVisoftLabels(filename)

strData = readmatrix(filename,'OutputType','string');
doubleData = str2double(strData(strData(:,3) == 'Alu'));
data.Alu = unique(doubleData);
data.Silence = str2double(strData((strData(:,3) == 'Silence'),1:2));

end

