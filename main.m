%% MEDICAL IMAGES ASSIGNMENT - TOPIC 1 - FBP ALGORITHM
% Loading the phantom
phant = phantom();

% Filters: None, Ram-Lak, Shepp-Logan and Hann
filtertype = {'none','Ram-Lak','Shepp-Logan','Hann'};

% Projection angles: from 0 to 180° with 20, 40, 70 or 180 values in between
thetaspan = [20,40,70,180];

% Noise: Gaussian with zero-mean and variance 1e-2, 1e-3, 1e-4 or no noise
noisevar = [1e-2, 1e-3, 1e-4, 0];

% Initializing cells where the projections and reconstructions are gonna be
% stored. We also keep track of the mean squared error 
radons = cell(4,4,4);
reconst = cell(4,4,4);
mse = zeros(4,4,4);

for f = 1:4 % Index 'f' identifies the Filters
    for t = 1:4 % Index 't' identifies the number of Theta projections
        for n = 1:4 % Index 'n' identifies the amount of Noise
            theta = linspace(0,180, thetaspan(t));
            % Noise addition
            noisyphant = imnoise(phant, 'gaussian',0, noisevar(n));
            % Radon transform, saved in the cell array
            radons{f,t,n} = radon(noisyphant, theta);
            
            % Reconstruction with 'iradon'
            reconst{f,t,n} = iradon(radons{f,t,n}, theta, 'linear', filtertype{f},1, size(phant,1));
            
            % MSE calculation
            mse(f,t,n) = mean((reconst{f,t,n}-noisyphant).^2, 'all');
        end
    end
end
% Identifying the best reconstruction as the one with the least MSE
[bestf, bestt, bestn] = ind2sub([4 4 4], find(mse == min(mse, [], 'all')));
% Identifying the worst reconstruction as the one with the highest MSE
[worstf, worstt, worstn] = ind2sub([4 4 4], find(mse == max(mse, [], 'all')));

%% VISUALIZATION OF THE BEST RECONSTRUCTED IMAGE
figure('Name','Best reconstruction');
subplot(121); imshow(phant, [],'InitialMagnification','fit');
title('Original Image');
subplot(122); imshow(reconst{bestf, bestt, bestn}, [],'InitialMagnification','fit');
title(['Best Reconstructed Image: Filt = ' filtertype{bestf} '  Proj = ' num2str(thetaspan(t))]);

%% SINOGRAM VISUALIZATION
% UNCOMMENT IF YOU WANT TO SEE THE SINOGRAMS
figure('Name','Radon Transforms');
for t=1:4
    for n=1:4
        subplot(4,4,t+4*n-4); imagesc(radons{1,t,n}); colormap parula
        title(['TH = ' num2str(thetaspan(t)) 'NV = ' num2str(noisevar(n))]);
    end
end

%% IRADON RECONSTRUCTION VISUALIZATION
% UNCOMMENT IF YOU WANT TO SEE THE RECONSTRUCTIONS FROM IRADON
for f = 1:4
    figure('Name',['Reconstruction - Filter: ' filtertype{f}]);
    for t=1:4
        for n=1:4
            subplot(4,4,t+4*n-4); imagesc(reconst{f,t,n}); colormap gray
            title(['TH = ' num2str(thetaspan(t)) 'NV = ' num2str(noisevar(n))]);
        end
    end
end

%% RECONSTRUCTON WITH THE IMPLEMENTED FBP ALGORITHM
rec_best = fbp(radons{bestf, bestt, bestn}, linspace(0,180,thetaspan(bestt)));
% The implemented function 'fbp' is at the bottom of the script

% Setting all the pixels in the background to 0 eliminate the star
% artifacts + CLEAN and SPUR operation to further reduce noise
rec_best(bwmorph(~imbinarize(rec_best),'clean','spur')) = 0; 

% Selecting the ROI in the original phantom
phant_bin = imbinarize(phant);
bb = regionprops(phant_bin,'BoundingBox'); bb = floor(bb.BoundingBox);
roi_phant = phant(bb(2):(bb(2)+bb(4)),bb(1):(bb(1)+bb(3)));

% Selecting the ROI in the iradon reconstruction
iradon_bin = imbinarize(reconst{bestf, bestt, bestn});
bb = regionprops(iradon_bin,'BoundingBox'); bb = floor(bb.BoundingBox);
roi_iradon = reconst{bestf, bestt, bestn}(bb(2):(bb(2)+bb(4)),bb(1):(bb(1)+bb(3)));

% Selecting the ROI in the implemented FBP reconstruction
rec_bin = imbinarize(rec_best); 
bb = regionprops(rec_best,'BoundingBox'); bb = floor(bb(1).BoundingBox);
roi_rec = rec_best(bb(2):(bb(2)+bb(4)),bb(1):(bb(1)+bb(3)));

% Setting all the pixels in the background to 0 eliminate the star
% artifacts
roi_rec(~imbinarize(roi_rec)) = 0;

% Median filtering to remove background noise
roi_rec_unfilt = roi_rec;
roi_rec = medfilt2(roi_rec,[3 3]);

% To calculate the MSE, the images must have the same size (there is
% probably 1/2 pixel difference in the sizes)
[sizer, sizec] = size(roi_phant);
roi_iradon = imresize(roi_iradon,[sizer, sizec]);
roi_rec = imresize(roi_rec,[sizer, sizec]);
roi_rec_unfilt = imresize(roi_rec_unfilt,[sizer, sizec]);

% Calculating the MSE in the 2 cases
mse_iradon = mean((roi_phant-roi_iradon).^2,'all');
mse_fbp = mean((roi_phant-roi_rec).^2,'all');
mse_fbp_iradon = mean((roi_iradon-roi_rec).^2,'all');

figure('Name','Reconstruction comparison');
subplot(131); imagesc(roi_phant); title('Original image'); 
subplot(132); imagesc(roi_iradon); title('Reconstructed with iradon');
subplot(133); imagesc(roi_rec); title('Reconstructed with FBP');
colormap gray

%% EXTENSION TO 3D PHANTOM
load('phant3d.mat');
theta = 0:180; % DeltaTHETA = 1°
slices = size(phant3d,3); % Number of slices 
radon3d = []; % Initializations
reconst3d = [];

% The radon transform is applied to each slice (this step takes a few
% seconds)
for s = 1:slices
    radon3d = cat(3, radon3d, radon(phant3d(:,:,s),theta));
    reconst3d = cat(3, reconst3d, fbp(radon3d(:,:,s),theta));
end

figure('Name','Radon transform - 3D');
montage(mat2gray(radon3d),'BorderSize',1,'BackgroundColor','white');
title('Radon transform of each slice of the 3D phantom');

figure('Name','Original 3D phantom'); montage(mat2gray(phant3d));
title('Original 3D phantom');

% ROI selection and background removal
reconst3droi = reconst3d(50:200,70:180,:);
reconst3droi(~imbinarize(reconst3droi)) = 0;

figure('Name','FBP-reconstructed 3D phantom'); montage(mat2gray(reconst3droi));
title('FBP-Reconstructed 3D phantom');

%% OUR FILTERED BACKPROJECTION ALGORITHM
function rec = fbp(sinog, theta)
%FBP(sinog, theta) applies the Filtered-Back-Projection algorithm to
%   reconstruct an image form the sinogram 'sinog', acquired at the
%   projection angles specified in the vector 'theta'

% Grayscale conversion
sinog = mat2gray(sinog); 

% Defining Ram-Lak filter
f = linspace(-1, 1, size(sinog,1))';
filter = fftshift(abs( f ));

% Filter is extended to cover all the projections
filter_mat = repmat(filter, [1 length(theta)]);
% FFT of the sinogram
sinog_freq = fft(sinog);
% Filtering in the frequency domain
sinog_filtered_f = sinog_freq.*filter_mat;
% Coming back to the space domain
sinog = real(ifft(sinog_filtered_f));
% The size of the image is the maximum one considering the inclination
% of the projection
N = size(sinog, 1);
Nmax = size(imrotate(zeros(N),45),1);
rec = zeros(Nmax,Nmax,length(theta));

for i = 1:size(sinog,2)
    % "Extrusion" of the projection
    current = repmat(sinog(:,i), [1,N]);
    % Rotation of angle theta (+90° so that the final image is aligned
    % with the original one);
    current_rotated = imrotate(current, theta(i)+90);
    % Padding with zeros
    current_padded = padarray(current_rotated,[floor((Nmax-size(current_rotated,1))/2), ...
        floor((Nmax-size(current_rotated,2))/2)]);
    if size(current_padded,1)<Nmax
        current_padded = padarray(current_padded,[1 1],'post');
    end
    rec(:,:,i) = current_padded;
end
rec = sum(rec,3);
end

