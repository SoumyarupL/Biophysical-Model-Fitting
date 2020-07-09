function [pixel] = fitting(fmel, fblood, model, T_RAW2XYZ, x)
    Sr = x(1,:);
    Sg = x(2,:);
    Sb = x(3,:);
    e = x(4,:);
    e = e./sum(e(:));
    %-------------------------- SKIN MODEL --------------------------
    [R_total] = skinModel(fmel,fblood,model);
    % Image Formation
    skincolour = R_total.*(e.*3); 
    rCh = sum(skincolour.*Sr);
    gCh = sum(skincolour.*Sg);
    bCh = sum(skincolour.*Sb);
    % raw image
    Iraw = cat(3, rCh,gCh,bCh);
    % Colour Transformation pipeline for visiualisation (see our BMVC 2019 paper)
    %1. White Balance
    [lightcolour] = computelightcolour(e,Sr,Sg,Sb); 
    [ImwhiteBalanced] = WhiteBalance(Iraw,lightcolour);
    % 2.from raw To RGB 
    [ pixel ] = reshape(fromRawTosRGB(ImwhiteBalanced,T_RAW2XYZ), [1 3]);
end