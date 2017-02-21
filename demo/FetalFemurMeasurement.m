function [starPoint,endPoint,SKL_ROI,cropImg]=FetalFemurMeasurement(USImg,USInfo,fanMask)
%This function measure fetal femur length for obstrict ultrasound image
%  input USimg-image data
%        USInfo- dicominfo
%        Mask- remove backgournd mask
% output startPoint endPoint-the end point of fetal femur
% create by Frank Zhao 05/13/2016

%pre_processing start
    x0 = USInfo.SequenceOfUltrasoundRegions.Item_1.RegionLocationMinX0;
    y0 = USInfo.SequenceOfUltrasoundRegions.Item_1.RegionLocationMinY0;
    x1 = USInfo.SequenceOfUltrasoundRegions.Item_1.RegionLocationMaxX1;
    y1 = USInfo.SequenceOfUltrasoundRegions.Item_1.RegionLocationMaxY1;
    figure,subplot(121);,imshow(USImg);
USImg=USImg(y0+1:y1,x0+1:x1,:);
[Height,Width,Slice]=size(USImg);
subplot(122);imshow(USImg);
 if Slice==3   
    IRed=USImg(:,:,1);
    IGreen=USImg(:,:,2);
    IBlue=USImg(:,:,3);
    IGray=IRed;
    IGrayLabel=bwlabel(IGray);
    CentValue=IGrayLabel(round(Height/2),round(Width/2));
    figure,subplot(221);imshow(IGrayLabel,[]);title('labeled');
    IGray(IGrayLabel~=CentValue)=0;
    IGrayLabel(IGrayLabel~=CentValue)=0;
    IGrayLabel=logical(IGrayLabel);
     subplot(222);imshow(IGrayLabel,[]);title('convert to binary');
    IGrayLabel=imfill(IGrayLabel,'holes');  
      subplot(223);imshow(IGrayLabel,[]); title('fill the holes');
    %IRedLabel Fun area;    
    Funedge=edge(IGrayLabel,'canny');
        subplot(224);imshow(Funedge,[]);title('canny');

    
    [rfan,cfan]=find(IGrayLabel>0);
    middleRow=round((rfan(1)+rfan(end))/2);
    middleCol=round((cfan(1)+cfan(end))/2);
    l=find(IGrayLabel(:,middleCol)>0);
    upperBoundary=round(l(1)+(l(end)-l(1))*0.13);
    bottomBoundary=round(l(1)+(l(end)-l(1))*0.75);

         Funedge(upperBoundary,:)=1;
    Funedge(bottomBoundary,:)=1;
    
%       figure,imshow(Funedge,[]);title('Image upper and lower boundary (0.13~0.75)') ;
     
    Mask=false(size(IGray));
    Mask(IRed==255 & IBlue==0)=1;
        

    
    if sum(Mask(:))==0
        RemoveIGray=IGray;
    else
        se=strel('disk',1);
        Mask1=imdilate(Mask,se);
        RemoveIGray=uint8(FastInpaint(IGray,Mask1,35));        
    end
 else
   IGrayLabel=bwlabel(USImg);
  figure,imshow(USImg,[]);
   IGray=USImg;
   figure,imshow(IGrayLabel,[]);
    CentValue=IGrayLabel(round(Height/2),round(Width/2));
    IGray(IGrayLabel~=CentValue)=0;
    IGrayLabel(IGrayLabel~=CentValue)=0;
    IGrayLabel=logical(IGrayLabel);
    IGrayLabel=imfill(IGrayLabel,'holes');
    %IRedLabel Fun area;
    Funedge=edge(IGrayLabel,'canny');
    [rfan,cfan]=find(IGrayLabel>0);
    middleRow=round((rfan(1)+rfan(end))/2);
    middleCol=round((cfan(1)+cfan(end))/2);
    l=find(IGrayLabel(:,middleCol)>0);
    upperBoundary=round(l(1)+(l(end)-l(1))*0.13);
    bottomBoundary=round(l(1)+(l(end)-l(1))*0.75);
%      Funedge(upperBoundary,:)=1;
    Funedge(bottomBoundary,:)=1;
    RemoveIGray=IGray;
 end
 
 %pre-processing end
 
cropImg=RemoveIGray;
%  literation OTSU 
 IBinROI=RemoveIGray;
    
H=fspecial('gaussian',5,2);
blurI=imfilter(RemoveIGray,H,'replicate');
    figure,subplot(121),imshow(IBinROI,[]);title('IBinROI');
       subplot(122),imshow(blurI,[]);title('Gaussian (5,2) Blur');
    %select the ROI region
M=blurI(IGrayLabel>0);
T2=graythresh(M)*255;

b=IBinROI;
b(b<T2)=0;
b(b>=T2)=255;
b=logical(b);   
    figure,subplot(221),imshow(b,[]);title('Fisrt OTSU');

M1=M(M>T2);
T=graythresh(M1)*255;

a=IBinROI;
a(a<T)=0;
a(a>=T)=255;
a=logical(a);   
    subplot(222),imshow(a,[]);title('Second OTSU');

T=(T/2+T2/2);
M2=M1(M1>T);
AdpThreshold=graythresh(M2)*255;
 
IBinROI(IBinROI<AdpThreshold)=0;
IBinROI(IBinROI>=AdpThreshold)=255;
IBinROI=logical(IBinROI);   
    subplot(223),imshow(IBinROI,[]);title('Third OTSU');

IBinROI=bwareaopen(IBinROI,270);   
     subplot(224),imshow(IBinROI,[]);title('Remove small object (pixels<270)');

%sine cut the topboundary
for i=1:Width
    A=round(sin(pi*(mod(i,36)/18))*10);
    IBinROI(1:A+upperBoundary+10,i)=0;
end
figure,subplot(221),imshow(IBinROI,[]);title('Cut top');

%remove the boundary block
se=strel('disk',1);
Funedge=imdilate(Funedge,se);
[LBC,N]=bwlabel(IBinROI);

a=IBinROI+Funedge;

subplot(222),imshow(a,[]);title('Cut top & Remove boundary block');

LBCAndFunedge=LBC.*Funedge;
boundaryNumber=unique(LBCAndFunedge(:));
for i=2:length(boundaryNumber)
   LBC(LBC==boundaryNumber(i))=0;   
end


% C is binary mask for high confidence region
HighConfROI=logical(LBC);



SeHor=strel('line',7,0);
HighConfROI=imclose(HighConfROI,SeHor);

subplot(223),imshow(HighConfROI,[]);title('High Confidence ROI after imclose using line(7,0)');

[labelC,CN]=bwlabel(HighConfROI);

subplot(224),imshow(labelC,[]);title('Labeled HighConfROI');

%extra skeleton 20 control the power of skrink
 thinedImage = bwmorph(HighConfROI,'skel',inf); 
 thinedImageshrink=bwmorph(thinedImage,'spur',20);
 

 figure,imshow(thinedImage,[]);title('image skeleton');
 figure,imshow(thinedImageshrink,[]);title('image spur');
 
 SKL_ROI=bwareaopen(thinedImageshrink,20);
  figure,imshow(SKL_ROI,[]);title('open using 20 after spur');
  
 if sum(SKL_ROI(:))==0
     thinedImageshrink=bwmorph(thinedImage,'spur',10);
     SKL_ROI=bwareaopen(thinedImageshrink,10);
      if sum(SKL_ROI(:))==0
          SKL_ROI=thinedImage;
      end
 end
 
 %debug =1
 [leftPoint, rightPoint,linestruct ] = SkeletonAnalysis( SKL_ROI,thinedImage,1);
 
 %end point analysis
 if leftPoint>0
%      if linestruct.theta~=0 && linestruct.rho~=0
         X=[leftPoint(2); rightPoint(2)];
         Y=[leftPoint(1); rightPoint(1)];
         % fit the reference line, 用fit主要是代码量偷懒
         lineFunction=fit(X,Y,'poly1');
         %drawLineimage用于画reference line
         drawlineimage=false(size(blurI));
         for i=1:size(blurI,2)
                drawlineimage(round(lineFunction.p1*i+lineFunction.p2),i)=1;
         end    
         

         
         %1 在指向线一定范围（上下高度）内进行连接
        %2 不能和标定段在X方向上有重叠 （先获取标定段）
        [labelC,CN]=bwlabel(HighConfROI);
        
%                  figure,imshow(drawlineimage+labelC,[]),title('drawline using poly1');
        %trustOne is the leftpoint locats in
        tagetnum=labelC(leftPoint(1),leftPoint(2));
        trustOne=labelC;
%                          figure,imshow(trustOne,[]),title('trustOne1');
        trustOne(trustOne~=tagetnum)=0;
        trustOne=logical(trustOne);
        [~,cTrust]=find(trustOne);
%                          figure,imshow(trustOne,[]),title('trustOne2');
        for i=1:CN
            [~,cTemp]=find(labelC==i);
         % newlines(i).point1(1)<centerRange(end) && newlines(i).point2(1)>centerRange(1) 
         if i~=tagetnum;
             if cTemp(1)<cTrust(end) && cTemp(end)>cTrust(1)
                labelC(labelC==i)=0;
             end
         end
        end
        lineTopWhite=drawlineimage;
        for i=1:size(blurI,2)
            A=drawlineimage(:,i);
            idxW=find(A>0);
            A(1:idxW(end))=1;
            lineTopWhite(:,i)=A;
        end
        
        lineTopWhite=logical(lineTopWhite);
%          figure,imshow(lineTopWhite),title('lineTopWhite');
        
        shrikROI=trustOne;
        shrikROI(lineTopWhite)=0;
         strLength=4;
         
%         figure,imshow(shrikROI),title('shrikROI');
        %use imdilate to find the ROI
        while sum(shrikROI(:))>0
        strLength=strLength+2;
        %line 90 degree strel to find the ROI
        seA=strel('line',strLength,90);
        lineTopWhite1=imdilate(lineTopWhite,seA);
        shrikROI=trustOne;
        shrikROI(lineTopWhite1)=0;
        end
        
%         figure,imshow(lineTopWhite1),title('lineTopWhite1');
        
        se=strel('line',strLength+1,90);
        lineTopWhite1=imdilate(lineTopWhite,se);
        labelC1=labelC;
        labelC1(lineTopWhite1)=0;

        
        %not meet this
        idx_lowThanTrust=find(labelC1>0);
        if ~isempty(idx_lowThanTrust)
            lowThanTrust=labelC(idx_lowThanTrust);
                    lowThanTrust=unique(lowThanTrust);
                    for i=1:length(lowThanTrust)
                        labelC(labelC==i)=0;
                    end 
                   HighConfROI=logical(labelC)+trustOne;
        else
            HighConfROI=trustOne+logical(labelC);
            
        end

         %在HighConfROI上操作，
         %1 imclose 连接
         %2 去掉不连接的
         %3 在斜率方向上的长宽比，过大，使用直线通过上的点。正常，是用两个端点
         
         %get a centblock as operator for imclose
         Operator_radius=13;
         centerline=round(size(blurI,2)/2);
         WhiteIdx=find(drawlineimage(:,centerline)==1);
         se=drawlineimage(WhiteIdx-Operator_radius:WhiteIdx+Operator_radius,centerline-Operator_radius:centerline+Operator_radius);
         HighConfROI_connect=imclose(HighConfROI,se);
        
 end
%  HighConfROI_connect
    TureROI=bwlabel(HighConfROI_connect);
    Femurskel=TureROI(leftPoint(1),leftPoint(2));  
    TureROI(TureROI~=Femurskel)=0;
    
%      figure,imshow(TureROI),title('TureROI');
    
    %normal direction -1/lineFunction.p1
%     k=-1/lineFunction.p1;
%     TestA=false(size(blurI));
%     for i=1:size(blurI,2)
%           TestA(round(k*i+lineFunction.p2),i)=255;
%     end


   k=-1/lineFunction.p1; %vertical direction
   Stand_Aspect_Ratio=11;%8;
BoneDirection=false(size(blurI));

for i=1:size(blurI,2)
    BoneDirection(round(lineFunction.p1*i+lineFunction.p2),i)=1;
end

figure,imshow(BoneDirection+TureROI),title('BoneDirection1');

[~,c]=find(TureROI>0);
BoneDirection(:,1:c(1)-1)=0;
BoneDirection(:,c(end)+1:end)=0;

[r,c]=find(BoneDirection>0);

%get the segment
figure,imshow(BoneDirection),title('BoneDirection2');
% beg=[r(1) c(1)];
% Ed=[r(end) c(end)];
% instcept_beg=r(1)-k*c(1);
% instcept_End=r(end)-k*c(end);
[Nr,idx,~] =unique(r,'stable');
Nc=c(idx);
isLarge=0;
smallest=200;
for j=1:length(c)
    NormalBoneDirection=false(size(blurI));
    intercept=r(j)-k*c(j);
    for i=1:size(blurI,1)
        col=round((i-intercept)/k);
        if col<1
            continue;
        end
        if col < size(blurI,2)
            NormalBoneDirection(i,col)=1;
        else
            break;
        end
    end
    NormalBoneDirection=logical(NormalBoneDirection) .* logical(TureROI);
    
%     figure,imshow(NormalBoneDirection),title('NormalBoneDirection1');
    
    if sum(NormalBoneDirection(:))>1
        [RN,CN]=find(NormalBoneDirection);
        Aspect_Ratio=sqrt((c(end)-c(1))^2+(r(end)-r(1))^2)/sqrt((RN(end)-RN(1))^2+(CN(end)-CN(1))^2);
        if Aspect_Ratio<smallest
            smallest=Aspect_Ratio;
        end
        if Aspect_Ratio>Stand_Aspect_Ratio
            isLarge=0;
        else
            isLarge=1;
            break;
        end
    end
end  

   
   se1=strel('disk',5);
   if isLarge
      BoneDirection=false(size(blurI));
      for i=1:size(blurI,2)
            BoneDirection(round(lineFunction.p1*i+lineFunction.p2),i)=255;
      end
       BoneDirection=imdilate(BoneDirection,se1);
       TureROI=BoneDirection&TureROI;
   end
        figure,imshow(BoneDirection),title('BoneDirection result');
         figure,imshow(TureROI),title('TureROI');
[r,c]=find(TureROI>0);
starPoint=[r(1),c(1)];
endPoint=[r(end),c(end)];
showresult=blurI;
showresult(r(1)-2:r(1)+2,c(1)-2:c(1)+2)=255;
showresult(r(end)-2:r(end)+2,c(end)-2:c(end)+2)=255;

 end
%%
function [leftPoint, rightPoint,linestruct ] = SkeletonAnalysis( SKLImg,thinedImage,debug)
[ThinSkelLabel,ThinSkelNum]=bwlabel(SKLImg);
figure,imshow(ThinSkelLabel);title('ThinSkelLabel')
% if ThinSkelNum>1
    [H,T,R] = hough(SKLImg,'RhoResolution',3,'Theta',-90:1:89);
    H(:,21:160)=0; % +-20degree only
   
    %debug output
  if debug>0
    figure,imshow(imadjust(mat2gray(H)),'XData',T,'YData',R,'InitialMagnification','fit');
    xlabel('\theta angle between x-axis'), ylabel('\rho');
    axis on, axis normal, hold on;
    colormap(hot);
  end
   
    P  = houghpeaks(H,3);
    %---debug output
  if debug>0
    figure,imshow(H,[],'XData',T,'YData',R,'InitialMagnification','fit');
    xlabel('\theta'), ylabel('\rho');
    axis on, axis normal, hold on;
    plot(T(P(:,2)),R(P(:,1)),'s','color','white');
  end  
  
    lines = houghlines(SKLImg,T,R,P,'FillGap',5,'MinLength',20);

    if isempty(lines)
        P  = houghpeaks(H,4);
        lines = houghlines(SKLImg,T,R,P,'FillGap',10,'MinLength',15);
    end
    %line structure analysis
    if debug>0
        figure,imshow(SKLImg,[]),hold on
        for k=1:length(lines)
            xy=[lines(k).point1; lines(k).point2];
            plot(xy(:,1),xy(:,2),'LineWidth',2,'Color',[1 0 0]);

        end
    end

%merge lines theta diff<5, rho diff<30
Theta_diff=5;
Rho_diff=30;
GroupMark=zeros(length(lines),1);
% GroupMark group index for all lines, e.g. (1 2 2 1) means 4 lines divide
% into 2 groups. line 1 and 4 belong to gourp 1, line 2 and 3 belong to
% group 2.
GroupMark(1,1)=1;
for i=2:length(lines)
    currenctGroupNum=max(GroupMark);
    for j=1:i-1
        if abs(lines(i).theta-lines(j).theta)< Theta_diff && abs(lines(i).rho-lines(j).rho)< Rho_diff
            GroupMark(i,1)=GroupMark(j,1);
            break;
        end
        if j==i-1
        GroupMark(i,1)=currenctGroupNum+1;
        end
    end    
end
Number_newlines=length(unique(GroupMark));
newlines=lines;
newlines(Number_newlines+1:end)=[];

    if debug>0
        figure,imshow(SKLImg,[]),title('new lines'),hold on
        for k=1:length(newlines)
            xy=[newlines(k).point1; newlines(k).point2];
            plot(xy(:,1),xy(:,2),'LineWidth',2,'Color',[1 0 0]);

        end
    end

%longest line on the different groups
for i=1:Number_newlines
    Groupidx=find(GroupMark==i);
    LineLength=0;
    longest_Idx=1;
    for j=1:length(Groupidx)
        xys=double([lines(Groupidx(j)).point1; lines(Groupidx(j)).point2]);
        currrentlength=sqrt((xys(1,1)-xys(2,1))^2+(xys(1,2)-xys(2,2))^2);
        if LineLength<currrentlength
            LineLength=currrentlength;
            longest_Idx=Groupidx(j);
        end
    end
    newlines(i)=lines(longest_Idx);
end

% select the buttom one of new lines. newlines(1).point1 [421   186] 
% 421 is x coordinate(horizontal) 186 is Y (vertical)
% x_length: [newlines(1).point1(1):newlines(1).point2(1)]
%find the lines have same coorindate(overlap) on horizontal direction
    if length(newlines)==1
        leftPoint=[newlines(1).point1(2) newlines(1).point1(1)]; % [r c]
        rightPoint=[newlines(1).point2(2) newlines(1).point2(1)]; % [r c]
        linestruct=newlines(1);
    else
        %find the x range on the center
        %define a cent mask
        HorCenter=round(size(SKLImg,2)/2);
        centerRadius=75;
        centerRange=HorCenter-centerRadius:HorCenter+centerRadius;
        largestVerCent=0;
        seletionMark=0;
        for i=1:length(newlines)
            %have overlap? A<End && B>Start
            if newlines(i).point1(1)<centerRange(end) && newlines(i).point2(1)>centerRange(1)
                %VerCent Largest
                VerCent=(newlines(i).point1(2)+newlines(i).point2(2))/2;
                if VerCent>largestVerCent
                    largestVerCent=VerCent;
                    seletionMark=i;
                end

            end
        end
        if 0==seletionMark
            leftPoint=-1;
            rightPoint=-1;
            linestruct=-1;
            return;
        end
        leftPoint=[newlines(seletionMark).point1(2) newlines(seletionMark).point1(1)]; % [r c]
        rightPoint=[newlines(seletionMark).point2(2) newlines(seletionMark).point2(1)]; % [r c]
        linestruct=newlines(seletionMark);
    end

% else
%         [Skelr,Skelc]=find(SKLImg>0); 
%        leftPoint=[Skelr(1),Skelc(1)];
%        rightPoint=[Skelr(end),Skelc(end)];  
%        linestruct.point1=[leftPoint(2) leftPoint(1)];
%        linestruct.point2=[rightPoint(2) rightPoint(1)];
%        linestruct.theta=0;
%        linestruct.rho=0;
% end

end
 %%
function k=FastInpaint(I,mask,liter)%(I,mask,liter,bar)
%Image inpainting by Manuel M. Oliveira's method Fast Digital Image
%Inpainting 
% 09-10-2007 ZhaoMing

%parameters:
%     I------image to be inpainted
%  mask------the noise mask which is a binary image
% liter------literation times
%   bar______diffusion barrier
% I=im2double(I);
I=double(I);
diffker=[0.073235 0.176765 0.0732325; 0.17675 0 0.17675; 0.073235 0.176775 0.073235]
%There is two kinds of diffusion kernels in Olivera's article I use the
%first one[0.073235 0.176765 0.0732325; 0.17675 0 0.17675; 0.073235 0.176775 0.073235]. 
%Another is [0.125 0.125 0.125; 0.125 0 0.125; 0.125 0.125 0.125].
if (size(I, 3) == 3)  %Color image process
    r=FInpaintGray(I(:,:,1),mask,liter,diffker);
    g=FInpaintGray(I(:,:,2),mask,liter,diffker);
    b=FInpaintGray(I(:,:,3),mask,liter,diffker);
    k=cat(3,r,g,b);
else                  %Gray image process
    k=FInpaintGray(I,mask,liter,diffker);
end
end
%-----------------------------------------------------------------------%
function g=FInpaintGray(I,mask,liter,diffker)
    Im=I;
    logmask=im2bw(mask);
    revmask=not(logmask);
    if nargin == 3
        for i=1:liter
           Im=imfilter(Im,diffker,'conv','replicate');
           Im=(Im.*logmask)+(I.*revmask);
         end
    else
        for i=1:liter
           Im=imfilter(Im,diffker,'conv','replicate');
           Im=(Im.*logmask)+(I.*revmask);
        end   
    end
    g=(Im.*logmask)+(I.*revmask);
    
end    




