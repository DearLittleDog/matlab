clear all;
close all;
clc;

fin=fopen('f:\processed.2048.2500.raw','r');
Z=fread(fin,2048*2500,'uint16');
Z=reshape(Z,2048,2500)';

for i=1:2500
    for j=1:2048
        B(i,j)=4096-double(Z(i,j));
    end
end
figure(1)
imshow(B,[]);

C=medfilt2(B,[3,3]);
