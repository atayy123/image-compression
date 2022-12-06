clear
%% apply FWT with Daubechies 8-tap filter to image harbour
% https://wavelets.pybytes.com/wavelet/db8/
db8 = [-0.00011747678400228192
0.0006754494059985568
-0.0003917403729959771
-0.00487035299301066
0.008746094047015655
0.013981027917015516
-0.04408825393106472
-0.01736930100202211
0.128747426620186
0.00047248457399797254
-0.2840155429624281
-0.015829105256023893
0.5853546836548691
0.6756307362980128
0.3128715909144659
0.05441584224308161
];

img = imread("images\harbour512x512.tif");
%imshow(img)
img = double(img(:));

y = fwt(img, db8);

i = 1:length(db8);
inv_db8 = db8(length(db8)-i+1);

result = inv_fwt(y, inv_db8);

img_reconst = reshape(result, [sqrt(length(result)), sqrt(length(result))]);
% figure
% imshow(uint8(img_reconst))

% compare the two images if they are equal, to test if our FWT
% implementation works.
diff = immse(result, img);
disp(diff)

%% scale 4 2D-FWT implementation
% To compute the wavelet coefficients at scale 4, we need 4 filter banks.

% filter bank 1
y1 = fwt(img, db8);

% get the low frequency part of the result
y1_low = y1(1:length(y1)/2);

% filter bank 2
y2 = fwt(y1_low, db8);
y2_low = y2(1:length(y2)/2);

% filter bank 3
y3 = fwt(y2_low, db8);
y3_low = y3(1:length(y3)/2);

% filter bank 4
y4 = fwt(y3_low, db8);

% get low and high frequency parts of the result
y4_low = y4(1:length(y4)/2);
y4_high = y4(length(y4)/2+1:end);

% plot the results
figure
imshow(uint8(reshape(y4_low, [sqrt(length(y4_low)), sqrt(length(y4_low))])))
title("Scale 4 Low frequency")

figure
imshow(uint8(reshape(y4_high, [sqrt(length(y4_high)), sqrt(length(y4_high))])))
title("Scale 4 High frequency")

%% Uniform Quantizer
stepsize = 2;

y4_low_quantized = stepsize * round(y4_low/stepsize);
y4_high_quantized = stepsize * round(y4_high/stepsize);

% plot the quantized results
figure
imshow(uint8(reshape(y4_low_quantized, [sqrt(length(y4_low)), sqrt(length(y4_low))])))
title("Scale 4 Low frequency quantized")

figure
imshow(uint8(reshape(y4_high_quantized, [sqrt(length(y4_high)), sqrt(length(y4_high))])))
title("Scale 4 High frequency quantized")

%% Distortion and Bit-Rate Estimation

% get images in the dataset
img_boat = double(imread("images\boats512x512.tif"));
img_peppers = double(imread("images\peppers512x512.tif"));
img_harbour = double(imread("images\harbour512x512.tif"));

% transform images to vectors
img_boat = img_boat(:);
img_peppers = img_peppers(:);
img_harbour = img_harbour(:);

% apply two step filterbank to the images
% filter bank 1
boat1 = fwt(img_boat, db8);
boat1_low = boat1(1:length(boat1)/2);
boat1_high = boat1(length(boat1)/2+1:end);

peppers1 = fwt(img_peppers, db8);
peppers1_low = peppers1(1:length(peppers1)/2);
peppers1_high = peppers1(length(peppers1)/2+1:end);

harbour1 = fwt(img_harbour, db8);
harbour1_low = harbour1(1:length(harbour1)/2);
harbour1_high = harbour1(length(harbour1)/2+1:end);

% filter bank 2
boat2 = fwt(boat1_low, db8);
boat2_low = boat2(1:length(boat2)/2);
boat2_high = boat2(length(boat2)/2+1:end);

peppers2 = fwt(peppers1_low, db8);
peppers2_low = peppers2(1:length(peppers2)/2);
peppers2_high = peppers2(length(peppers2)/2+1:end);

harbour2 = fwt(harbour1_low, db8);
harbour2_low = harbour2(1:length(harbour2)/2);
harbour2_high = harbour2(length(harbour2)/2+1:end);

d_b = zeros(1,10);
d_p = zeros(1,10);
d_h = zeros(1,10);

e_b = zeros(1,10);
e_p = zeros(1,10);
e_h = zeros(1,10);
for i = 0:9
    % apply quantization
    stepsize = 2^i;
    q_b2_low = stepsize * round(boat2_low/stepsize);
    q_b2_high = stepsize * round(boat2_high/stepsize);
    q_b1_high = stepsize * round(boat1_high/stepsize);

    q_p2_low = stepsize * round(peppers2_low/stepsize);
    q_p2_high = stepsize * round(peppers2_high/stepsize);
    q_p1_high = stepsize * round(peppers1_high/stepsize);

    q_h2_low = stepsize * round(harbour2_low/stepsize);
    q_h2_high = stepsize * round(harbour2_high/stepsize);
    q_h1_high = stepsize * round(harbour1_high/stepsize);

    % reconstruct images
    b1low = inv_fwt([q_b2_low; q_b2_high], inv_db8);
    b_recons = inv_fwt([b1low; q_b1_high], inv_db8);

    p1low = inv_fwt([q_p2_low; q_p2_high], inv_db8);
    p_recons = inv_fwt([p1low; q_p1_high], inv_db8);

    h1low = inv_fwt([q_h2_low; q_h2_high], inv_db8);
    h_recons = inv_fwt([h1low; q_h1_high], inv_db8);
    %figure
    %imshow(uint8(reshape(h_recons, [512, 512])))

    % calculate distortion
    d_b(i+1) = immse(b_recons, img_boat);
    d_p(i+1) = immse(p_recons, img_peppers);
    d_h(i+1) = immse(h_recons, img_harbour);

    % calculate the entropy (bit-rate estimation)
    e_b(i+1) = (4*entropy(q_b2_low) + 4*entropy(q_b2_high) + 2*entropy(q_b1_high));
    e_p(i+1) = (4*entropy(q_p2_low) + 4*entropy(q_p2_high) + 2*entropy(q_p1_high));
    e_h(i+1) = (4*entropy(q_h2_low) + 4*entropy(q_h2_high) + 2*entropy(q_h1_high));
end

% compare d and difference between wavelet coefficients and their quantized
% versions

q_dist = (4*immse(boat2_high, q_b2_high) + 2*immse(boat1_high, q_b1_high) + 4*immse(boat2_low, q_b2_low))/10;
disp(d_h(end))
disp(q_dist)

% The distortion between the original and reconstructed image is close to
% the average of the distortion of the wavelet coefficients and their
% quantized versions.

% calculate PSNR
p_b = 10*log10(255^2./d_b);
p_h = 10*log10(255^2./d_h);
p_p = 10*log10(255^2./d_p);

% plot the figures
figure
plot(d_b, e_b)
hold on
plot(d_h,e_h)
plot(d_p,e_p)
legend(["Boats", "Harbour", "Peppers"])
xlabel("Distortion")
ylabel("Rate")

figure
plot(e_b, p_b)
hold on
plot(e_h, p_h)
plot(e_p, p_p)
legend(["Boats", "Harbour", "Peppers"])
ylabel("PSNR")
xlabel("Rate")