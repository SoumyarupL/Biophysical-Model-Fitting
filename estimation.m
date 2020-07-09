addpath('matFiles');
load XYZspace.mat; % The CIE XYZ color space.
load rgbCMF; % databse for 28 cameras
% load the histological quantities and functional approximations.
load pheomelanin_ext;
load eumelanin_ext;
load deoxy_hemo_ext_coeff;
load oxy_hemo_ext_coeff;
S = [Sr;Sg;Sb];  %S : 3 x 33
T = (S'\XYZspace)'; 
T(1,:) = 0.3253.*T(1,:)./sum(T(1,:));
T(2,:) = 0.3425.*T(2,:)./sum(T(2,:));
T(3,:) = 0.3723.*T(3,:)./sum(T(3,:));
T_RAW2XYZ = reshape(T, [1 1 9]);
% save the histological as struct
[model] = preparedModel(pheomelanin_ext,eumelanin_ext,deoxy_hemo_ext_coeff,oxy_hemo_ext_coeff);
[Y] = CameraSensitivity(rgbCMF);
wavelength = 33;
Sr = reshape(double(Y(1:wavelength,1)),[1 33]);
Sg = reshape(double(Y(wavelength+1:wavelength*2,1)),[1 33]);
Sb = reshape(double(Y(wavelength*2+1:wavelength*3,1)),[1 33]);

load illumA; % spectra of illuminant A
e = reshape(double(illumA),[1 33]);
e = e./sum(e(:));

minmelanin = 0.013;  
maxmelanin = 0.43; 
minhemoglobin = 0.02; 
maxhemoglobin = 0.07;

xdata = [Sr;Sg;Sb;e];
J = imcrop(imread('image2.jpg'));
figure;
imshow(J);
fmel = zeros(size(J,1),size(J,2));
fblood = zeros(size(J,1),size(J,2));
[nRows, nCols] = size(J);
tic
for row = 1:nRows
    fm = zeros(1,size(J,2));
    fb = zeros(1,size(J,2));
    for col = 1:nCols/3
        ydata = reshape(double(J(row,col,:)), [1 3])./255;
        fun = @(x,xdata)fitting(x(1), x(2), model, T_RAW2XYZ, xdata);
        lb = [minmelanin, minhemoglobin];
        ub = [maxmelanin, maxhemoglobin];
        x0 = [0.2, 0.03];
        x = lsqcurvefit(fun, x0, xdata, ydata, lb, ub);
        fm(col) = x(1);
        fb(col) = x(2);
        fprintf("%f\n",x);
    end
    fmel(row,:) = fm;
    fblood(row,:) = fb;
end
toc

