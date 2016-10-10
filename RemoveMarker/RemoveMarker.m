function [ret_img] = RemoveMarker(dicom_info, img)
color_type = dicom_info.ColorType;

if strcmp(color_type, 'grayscale')
    ret_img = RemoveMarkerInGrayImg(dicom_info, img);
else
    if strcmp(dicom_info.Manufacturer, 'Carestream Health') && ...
            strcmp(dicom_info.ManufacturerModelName, 'Touch Ultra Sound System')
        ret_img = RemoveMarkerInCarestreamColorImg(dicom_info, img);
    else
        ret_img = RemoveMarkerInColorImg(dicom_info, img);
    end
end
%% convert to original color scheme
% ret_img = ConvertToOriginalColorScheme(dicom_info, ret_img);
%%
if false
    line_bound = 255 * ones(size(ret_img, 1), 1, size(ret_img, 3));
    img_region = GetImageRegion(dicom_info, img);
    ret_img = [img_region  line_bound ret_img];
    ret_img = insertText(ret_img, [1 1; size(img_region, 2)+ 1 1],{'before', 'after'});
end
ret_img = uint8(ret_img);
end

function [ret_img] = RemoveMarkerInGrayImg(dicom_info, img)
%threshold value 254, may need one table for different vendor and model, TBD
tv = 253;
%%
img_region = GetImageRegion(dicom_info, img);
mask = img_region > tv;
MN = [4 4];
SE = strel('rectangle',MN);
mask_d = imdilate(mask, SE);
e = (mask_d - mask) ~= 0;
mask((img_region == 1 | img_region == 0) & e) = 1;
ret_img= FillImageHole(img_region, mask, 'inpaint');
end

function [ret_img] = RemoveMarkerInRGBImg(dicom_info, img)

diff_level = 0.3;
img_region = GetImageRegion(dicom_info, img);
mask = CalMask(dicom_info, img, diff_level);
ret_img= FillImageHole(img_region, mask, 'inpaint');
end

function [ret_img] = RemoveMarkerInIndexedImg(dicom_info, img)
diff_level = 0.3;
img_region = GetImageRegion(dicom_info, img);
mask = CalMask(dicom_info, img, diff_level);
ret_img= FillImageHole(img_region, mask, 'inpaint');
end

function [ret_img] = RemoveMarkerInColorImg(dicom_info, img)
diff_level = 0.3;
img_region = GetImageRegion(dicom_info, img);
mask = CalMask(dicom_info, img, diff_level);
ret_img= FillImageHole(img_region, mask, 'inpaint');
end

function [ret_img] = RemoveMarkerInCarestreamColorImg(dicom_info, img)
%threshold value 254, may need one table for different vendor and model, TBD
tv = 253;
%%
img_region = GetImageRegion(dicom_info, img);
mask = img_region(:,:,1) > tv;
MN = [4 4];
SE = strel('rectangle',MN);
mask_d = imdilate(mask, SE);
e = (mask_d - mask) ~= 0;
% figure, imshow(e, [])
mask((img_region(:,:,1) == 1 | img_region(:,:,1) == 0) & e) = 1;
ret_img= FillImageHole(img_region, mask, 'inpaint');
end

function [img_region] = GetImageRegion(info, img)

x0 = info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMinX0;
y0 = info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMinY0;
x1 = info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMaxX1;
y1 = info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMaxY1;
img_region = img(y0+1:y1, x0+1:x1, :);
% img_region = img;
if strcmp(info.ColorType, 'indexed')
    %TBD, may need check info online
    n = 256;
    c = 3;
    if isfield(info, 'RedPaletteColorLookupTableDescriptor')
        des = info.RedPaletteColorLookupTableDescriptor;
    end
    n = des(1);
    map_ = zeros(n, c);
    map_(:,1) = single(bitand(info.RedPaletteColorLookupTableData, 255))/255.0;
    map_(:,2) = single(bitand(info.GreenPaletteColorLookupTableData, 255))/255.0;
    map_(:,3) = single(bitand(info.BluePaletteColorLookupTableData, 255))/255.0;
    img_region = ind2rgb(img_region, map_) * 255;
end

end

function [mask] = CalMask(info, img, diff_level)

img_region = GetImageRegion(info, img);
img_region_channels_mean = mean(img_region, 3);
max_diff = max(abs(single(img_region) - repmat(img_region_channels_mean, 1, 1, 3)), [], 3);
mask = max_diff > diff_level.*img_region_channels_mean;
MN = [4 4];
SE = strel('rectangle',MN);
mask_d = imdilate(mask, SE);
e = (mask_d - mask) ~= 0;
mask((img_region_channels_mean == 1 | img_region_channels_mean == 0) & e) = 1;
end

function [ret_img] = FillImageHole(img, mask, method)
z = size(img, 3);
ret_img = zeros(size(img));
if strcmp(method, 'regionfill')
    for i = 1:z
        ret_img(:,:,i) = regionfill(img(:,:,i), mask);
    end
elseif(strcmp(method, 'inpaint'))
    t = double(img);
    t(repmat(mask, 1, 1, z)) = NaN('double');
    for i = 1:z
        ret_img(:,:,i) = inpaint_nans(t(:,:,i), 3);
    end
end
end

function [ret_img] = ConvertToOriginalColorScheme(dicom_info, img)

if isfield(dicom_info, 'RedPaletteColorLookupTableDescriptor')
    map(:, 1) = single(bitand(dicom_info.RedPaletteColorLookupTableData, 255))/255.0;
    map(:, 2) = single(bitand(dicom_info.GreenPaletteColorLookupTableData, 255))/255.0;
    map(:, 3) = single(bitand(dicom_info.BluePaletteColorLookupTableData, 255))/255.0;
    ret_img = rgb2ind(img, map);
else
    ret_img = img;
end
end