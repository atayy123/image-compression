%% 8x8 DCT block

% calculate the A matrix
M = 8;
A = zeros(M,M);

for i = 1:M
    for k = 1:M
        if i == 1
            alpha = sqrt(1/M);
        else
            alpha = sqrt(2/M);
        end
        A(i,k) = alpha * cos(((2*(k-1)+1)*(i-1)*pi)/(2*M));
    end
end

% write functions for blockwise operation
dct = @(block_struct) A * block_struct.data * A';
inverse_dct = @(block_struct) A' * block_struct.data * A;

% print A values
disp(A)

%% apply DCT to an image
img = double(imread("images\airfield512x512.tif"));

dct_result = blockproc(img, [8 8], dct);


%% uniform quantizer, unlimited number of quantization steps
stepsize = 2;

quantized = stepsize * round(dct_result/stepsize); 

syms q(x)
q(x) = round(x/stepsize);
fplot(q)
grid on
xlabel("Input")
ylabel("Quantized output")

%% reconstruct the image and calculate the MSE

reconstructed = blockproc(quantized, [8 8], inverse_dct);
d = immse(reconstructed, img);
d_quantizer = immse(quantized, dct_result);

% distortion of the images and DCT coefficients are equal.

% figure
% imshow(uint8(reconstructed))

%% calculate rate and PSNR and plot the graph

% get images in the dataset
img_boat = double(imread("images\boats512x512.tif"));
img_peppers = double(imread("images\peppers512x512.tif"));
img_harbour = double(imread("images\harbour512x512.tif"));

% perform DCT to images
dct_boat = blockproc(img_boat, [8 8], dct);
dct_peppers = blockproc(img_peppers, [8 8], dct);
dct_harbour = blockproc(img_harbour, [8 8], dct);

% quantize the images with different step sizes
d_b = zeros(1,10);
d_p = zeros(1,10);
d_h = zeros(1,10);

e_b = zeros(1,10);
e_p = zeros(1,10);
e_h = zeros(1,10);

for i = 0:9
    stepsize =  2^i;

    quantized_b = stepsize * round(dct_boat/stepsize);
    quantized_p = stepsize * round(dct_peppers/stepsize);
    quantized_h = stepsize * round(dct_harbour/stepsize);

    % reconstruct the images
    b_rec = blockproc(quantized_b, [8 8], inverse_dct);
    p_rec = blockproc(quantized_p, [8 8], inverse_dct);
    h_rec = blockproc(quantized_h, [8 8], inverse_dct);
    
    % calculate the distortion
    d_b(i+1) = immse(quantized_b, dct_boat);
    d_p(i+1) = immse(quantized_p, dct_peppers);
    d_h(i+1) = immse(quantized_h, dct_harbour);

    % calculate the bit rate (in this case, the estimation is entropy,
    % because we assumed that we use the ideal code word length of VLC.)
    e_b(i+1) = entropy(uint8(b_rec));
    e_p(i+1) = entropy(uint8(p_rec));
    e_h(i+1) = entropy(uint8(h_rec));
%     e_b(i+1) = entropy(quantized_b);
%     e_p(i+1) = entropy(quantized_p);
%     e_h(i+1) = entropy(quantized_h);
end

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


