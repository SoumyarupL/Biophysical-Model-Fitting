function [mat] = reconstruct(fmel, fblood, x)
    addpath('matFiles');
    % load the histological quantities and functional approximations.
    load pheomelanin_ext;
    load eumelanin_ext;
    load deoxy_hemo_ext_coeff;
    load oxy_hemo_ext_coeff;
    load XYZspace.mat;
    [model] = preparedModel(pheomelanin_ext,eumelanin_ext,deoxy_hemo_ext_coeff,oxy_hemo_ext_coeff);
    Sr = x(1,:);
    Sg = x(2,:);
    Sb = x(3,:);
    e = x(4,:);
    e = e./sum(e(:));
    ar_row=size(fmel,1);
    ar_cl=size(fmel,2);
    SpectralReflectance = zeros(size(fmel,1),size(fmel,2),33);
     for row = 1:ar_row
         for col = 1:ar_cl
             m = fmel(row,col);
             h = fblood(row,col);
            [R_total] = skinModel(m,h,model); % compute spectral reflectance
             SpectralReflectance(row,col,:) =  R_total;
         end
     end 
    % Image Formation
    skincolour = SpectralReflectance.*reshape(e.*3,[1 1 33]); 
    rCh = sum(skincolour.*reshape(Sr,[1 1 33]),3);
    gCh = sum(skincolour.*reshape(Sg,[1 1 33]),3);
    bCh = sum(skincolour.*reshape(Sb,[1 1 33]),3);
    % raw image
    Iraw = cat(3, rCh,gCh,bCh);


    % Colour Transformation pipeline for visiualisation (see our BMVC 2019 paper)
    %1. White Balance
    [lightcolour] = computelightcolour(e,Sr,Sg,Sb); 
    [ImwhiteBalanced] = WhiteBalance(Iraw,lightcolour);
    % 2. Find T matrix from raw to xzy space
    S = [Sr;Sg;Sb];  %S : 3 x 33
    T = (S'\XYZspace)'; 
    T(1,:) = 0.3253.*T(1,:)./sum(T(1,:));
    T(2,:) = 0.3425.*T(2,:)./sum(T(2,:));
    T(3,:) = 0.3723.*T(3,:)./sum(T(3,:));
    T_RAW2XYZ = reshape(T, [1 1 9]);
    % 3.from raw To RGB 
    [ mat ] = fromRawTosRGB(ImwhiteBalanced,T_RAW2XYZ);
    imshow(mat);

end