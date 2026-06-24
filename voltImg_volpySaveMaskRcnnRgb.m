function voltImg_volpySaveMaskRcnnRgb(meanImgZ, corrImgZ, outBasePath)
%VOLTIMG_VOLPYSAVEMASKRCNNRGB Save 3-channel PNG for Mask R-CNN (VolPy / COCO-style input).
%   VolPy duplicates the mean in two channels and puts the correlation image in the third.
%   meanImgZ, corrImgZ: H x W (z-scored summary images).

lo = -3; hi = 3;
toU8 = @(z) uint8(255 * min(max((z - lo) / (hi - lo), 0), 1));
ch1 = toU8(meanImgZ);
ch2 = toU8(meanImgZ);
ch3 = toU8(corrImgZ);
rgb = cat(3, ch1, ch2, ch3);
outPng = [outBasePath, '_volpy_maskrcnn_rgb.png'];
imwrite(rgb, outPng);
fprintf('Saved Mask R-CNN style RGB for external inference: %s\n', outPng);
end
