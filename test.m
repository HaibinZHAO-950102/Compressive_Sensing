[a1, map]=imread('Kalman_modal_4.gif','frame','all') ;
[rows, columns, numColorChannels, numImages] = size(a1);
rgbImage = zeros(rows, columns, 3, numImages, 'uint8'); % Initialize dimensions.
for k = 1 : numImages
  thisFrame = a1(:,:,:, k);
  thisRGB = uint8(255 * ind2rgb(thisFrame, map));
  imshow(thisRGB);
  rgbImage(:,:,:,k) = thisRGB;
  drawnow;
end