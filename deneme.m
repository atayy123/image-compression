%% 8x8 DCT block
% calculate the A matrix
M = 8; % Block size
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
ent = @(block_struct) reshape(transpose(block_struct.data), 1, []); 
% prepare sorting the coefficients for entropy calculation

% print A values
disp(A)

%% apply DCT to an image
img = double(imread("images\airfield512x512.tif"));
dct_result = blockproc(img, [8 8], dct);

%% uniform quantizer, unlimited number of quantization steps
stepsize = 2;
quantized = stepsize * round(dct_result/stepsize); 

syms q(x)
q(x) = stepsize * round(x/stepsize);
figure
fplot(q)
grid on
xlabel("Input")
ylabel("Quantized output")

%% reconstruct the image and calculate the MSE
reconstructed = blockproc(quantized, [8 8], inverse_dct);
d = immse(reconstructed, img);
d_quantizer = immse(quantized, dct_result);
% distortion of the images and DCT coefficients are equal.

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
% distortion for different step sizes
d_b = zeros(1,10);
d_p = zeros(1,10);
d_h = zeros(1,10);
entro = zeros(1,10);

% calculate distortion and bite rate for step sizes 2^0 to 2^9
for i = 0:9
    stepsize =  2^i;

    quantized_b = stepsize * round(dct_boat/stepsize);
    quantized_p = stepsize * round(dct_peppers/stepsize);
    quantized_h = stepsize * round(dct_harbour/stepsize);

    % calculate the distortion as mse of the image
    d_b(i+1) = immse(quantized_b, dct_boat);
    d_p(i+1) = immse(quantized_p, dct_peppers);
    d_h(i+1) = immse(quantized_h, dct_harbour);

    [n,m] = size(quantized_h);
    sz = n*m;
    
    % transform coefficients of each block into a matrix of 64 rows
    ent_b = blockproc(quantized_b, [8 8], ent);
    ent_b = reshape(ent_b', 64, sz/64)';
    ent_p = blockproc(quantized_p, [8 8], ent);
    ent_p = reshape(ent_p', 64, sz/64)';
    ent_h = blockproc(quantized_h, [8 8], ent);
    ent_h = reshape(ent_h', 64, sz/64)';
    en_all = [ent_b; ent_p; ent_h];
    % sum over the entropies of all coefficient to get the total entropy of
    % a block
    for m = 1:64
        entro(i+1) = entro(i+1) + entropy(mat2gray(en_all(:,m)))/64;
    end
end

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