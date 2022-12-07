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

img = double(imread("images\harbour512x512.tif"));

[yll, yhl, ylh, yhh] = fwt2d(img, db8);

i = 1:length(db8);
inv_db8 = db8(length(db8)-i+1);

result = inv_fwt2d(yll, yhl, ylh, yhh, inv_db8);

% figure
% imshow(uint8(img_reconst))

% compare the two images if they are equal, to test if our FWT
% implementation works.
diff = immse(result, img);
disp(diff)

%% scale 4 2D-FWT implementation
% To compute the wavelet coefficients at scale 4, we need 4 filter banks.

% filter bank 1
[yll1, yhl1, ylh1, yhh1] = fwt2d(img, db8);

% filter bank 2 to the Low Low part of the FWT
[yll2, yhl2, ylh2, yhh2] = fwt2d(yll1, db8);

% filter bank 3
[yll3, yhl3, ylh3, yhh3] = fwt2d(yll2, db8);

% filter bank 4
[yll4, yhl4, ylh4, yhh4] = fwt2d(yll3, db8);

% plot the results
figure
tiledlayout(2,2)

nexttile
imshow(uint8(yll4))
title("LL")

nexttile
imshow(uint8(yhl4))
title("HL")

nexttile
imshow(uint8(ylh4))
title("LH")

nexttile
imshow(uint8(yhh4))
title("HH")

%% Uniform Quantizer
stepsize = 2;

yll4_q = stepsize * round(yll4/stepsize);
yhl4_q = stepsize * round(yhl4/stepsize);
ylh4_q = stepsize * round(ylh4/stepsize);
yhh4_q = stepsize * round(yhh4/stepsize);

% % plot the quantized results
% figure
% tiledlayout(2,2)
% title("Quantized Coefficients")
% 
% nexttile
% imshow(uint8(yll4_q))
% title("LL")
% 
% nexttile
% imshow(uint8(yhl4_q))
% title("HL")
% 
% nexttile
% imshow(uint8(ylh4_q))
% title("LH")
% 
% nexttile
% imshow(uint8(yhh4_q))
% title("HH")

%% Distortion and Bit-Rate Estimation

% get images in the dataset
img_boat = double(imread("images\boats512x512.tif"));
img_peppers = double(imread("images\peppers512x512.tif"));
img_harbour = double(imread("images\harbour512x512.tif"));

% apply 2d FWT to images
[b_ll, b_hl, b_lh, b_hh] = fwt2d(img_boat, db8);
[p_ll, p_hl, p_lh, p_hh] = fwt2d(img_peppers, db8);
[h_ll, h_hl, h_lh, h_hh] = fwt2d(img_harbour, db8);

d_b = zeros(1,10);
d_p = zeros(1,10);
d_h = zeros(1,10);

% e_b = zeros(1,10);
% e_p = zeros(1,10);
% e_h = zeros(1,10);
entro = zeros(1,10);
for i = 0:9
    % apply quantization
    stepsize = 2^i;
    q_b_ll = stepsize * round(b_ll/stepsize);
    q_b_hl = stepsize * round(b_hl/stepsize);
    q_b_lh = stepsize * round(b_lh/stepsize);
    q_b_hh = stepsize * round(b_hh/stepsize);

    q_p_ll = stepsize * round(p_ll/stepsize);
    q_p_hl = stepsize * round(p_hl/stepsize);
    q_p_lh = stepsize * round(p_lh/stepsize);
    q_p_hh = stepsize * round(p_hh/stepsize);

    q_h_ll = stepsize * round(h_ll/stepsize);
    q_h_hl = stepsize * round(h_hl/stepsize);
    q_h_lh = stepsize * round(h_lh/stepsize);
    q_h_hh = stepsize * round(h_hh/stepsize);

    % reconstruct images
    b_recons = inv_fwt2d(q_b_ll, q_b_hl, q_b_lh, q_b_hh, inv_db8);

    p_recons = inv_fwt2d(q_p_ll, q_p_hl, q_p_lh, q_p_hh, inv_db8);

    h_recons = inv_fwt2d(q_h_ll, q_h_hl, q_h_lh, q_h_hh, inv_db8);

    % calculate distortion
    d_b(i+1) = immse(b_recons, img_boat);
    d_p(i+1) = immse(p_recons, img_peppers);
    d_h(i+1) = immse(h_recons, img_harbour);

    % calculate the entropy (bit-rate estimation)
%     e_b(i+1) = (entropy(mat2gray(q_b_ll)) + entropy(mat2gray(q_b_hl)) + entropy(mat2gray(q_b_lh)) + entropy(mat2gray(q_b_hh)))/4;
%     e_p(i+1) = (entropy(mat2gray(q_p_ll)) + entropy(mat2gray(q_p_hl)) + entropy(mat2gray(q_p_lh)) + entropy(mat2gray(q_p_hh)))/4;
%     e_h(i+1) = (entropy(mat2gray(q_h_ll)) + entropy(mat2gray(q_h_hl)) + entropy(mat2gray(q_h_lh)) + entropy(mat2gray(q_h_hh)))/4;
    entro(i+1) = (entropy(mat2gray([q_b_ll q_p_ll q_h_ll])) + entropy(mat2gray([q_b_hl q_p_hl q_h_hl])) + entropy(mat2gray([q_b_lh q_p_lh q_h_lh])) + entropy(mat2gray([q_b_hh q_p_hh q_h_hh])))/4;
end

% compare d and difference between wavelet coefficients and their quantized
% versions

q_dist = (immse(b_ll, q_b_ll) + immse(b_hl, q_b_hl) + immse(b_lh, q_b_lh) + immse(b_hh, q_b_hh))/4;
disp(d_b(end))
disp(q_dist)

% The distortion between the original and reconstructed image is equal to
% the average of the distortion of the wavelet coefficients and their
% quantized versions.

% calculate PSNR
p_b = 10*log10(255^2./d_b);
p_h = 10*log10(255^2./d_h);
p_p = 10*log10(255^2./d_p);

% plot the figures
figure
plot(d_b, entro)
hold on
grid on
plot(d_h,entro)
plot(d_p,entro)
legend(["Boats", "Harbour", "Peppers"])
xlabel("Distortion")
ylabel("Rate")

figure
plot(entro, p_b)
hold on
grid on
plot(entro, p_h)
plot(entro, p_p)
legend(["Boats", "Harbour", "Peppers"])
ylabel("PSNR")
xlabel("Rate")