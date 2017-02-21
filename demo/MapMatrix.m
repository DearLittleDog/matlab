clear all;
close all;
clc;

smallpath='empty.raw';
smallimage='image.raw';

fin=fopen(smallpath,'r');
emptyimg=fread(fin,2520*3032,'uint16');
emptyimg=reshape(emptyimg,2520,3032)';

fin=fopen(smallimage,'r');
img=fread(fin,2520*3032,'uint16');
img=reshape(img,2520,3032)';

%original image
figure,imshow(emptyimg,[]);
% figure,imshow(img,[])

mapMatrix=ones(size(emptyimg));


xv= [458,1792,1775,440,458];
yv= [2396,2383,628,641,2396];
mapMask=poly2mask(xv,yv,3032,2520);
mapMask=mapMask.*2795;
mapMask=mapMask./emptyimg;
% mapMatrix=2795 /emptyimg;
figure,imshow(mapMask,[]);

result=img.*mapMask;
figure,imshow(result,[]);

resultid=fopen('result.raw','wb');
fwrite(resultid, result ,'uint16');

% 
% for x=1:3032
%     for y=1:2520
%         if inpolygon(x,y,xv,yv)
%             mapMatrix(x,y)=2795/emptyimg(x,y);
%         end
%     end
% end



