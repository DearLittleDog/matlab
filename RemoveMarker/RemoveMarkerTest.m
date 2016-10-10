clc
clear all
close all
folder = 'F:\Data\4 categories_alldatabases\URMCsave_Philips\';
id = strfind(folder, '\');
start = id(end - 1);
name = folder(start+1:end-1);
out_folder =strcat(folder, name,'_RemoveMarker_',date,'\');
mkdir(out_folder);
List = dir(folder);
toWrite = {'File','Marker', 'Success'};
for i = 1:size(List, 1)
    
    file = List(i).name;
    if isdir(strcat(folder, file)) || ~(isempty(strfind(file, 'dcm')) || isempty(strfind(file, 'DCM')))
        continue;
    end
    file_name = file;
    if ~isempty(strfind(file, '.'))
        file_name = file(1:end-4);
    end
    info = dicominfo(strcat(folder, file));
    img = dicomread(strcat(folder, file));
   
    ret_img = RemoveMarker(info, img);
    
    imwrite(ret_img, strcat(out_folder, file_name, '.bmp'));
%     figure, imshow(ret_img, []);
    toWrite = [toWrite; {file_name, 0, 1}];
end

xlswrite(strcat(out_folder,'Cases.xlsx'), toWrite);