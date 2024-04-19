function [ Ajoint, Abase, Fbase ] = bilateral( imflash, imambient,sigmag,ws,sigmab1)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    imambient = double(imambient);
    imflash = double(imflash);
    
%     sigmag = 5;
%     ws = 11;
%     sigmab1 = 50;
%     sigmab1 = (max(imambient(:)) - min(imambient(:)))*(255/10);



    gauss_mask = fspecial('gaussian',ws, sigmag);
    
    bias = floor(ws/2);
    flashpad = padarray(imflash,[bias bias],'replicate');
    ambientpad = padarray(imambient,[bias bias],'replicate');

    [h,w,c] = size(imflash);
    Ajoint = zeros(h,w,c);
    Abase = zeros(h,w,c);
    Fbase = zeros(h,w,c);
    for i = 1+bias:h-bias
        for j = 1+bias:w-bias
            amb_mask = ambientpad(i-bias:i+bias,j-bias:j+bias,:);
            flash_mask = flashpad(i-bias:i+bias,j-bias:j+bias,:);
            flash_diffmask = flash_mask-flashpad(i,j,:);
            amb_diffmask = amb_mask-ambientpad(i,j,:);
%             bil_mask_flash = exp(double(-1*((flash_diffmask/sigmab1).^2)/(2*(sigmab1^2))));      %/(sqrt(2*pi)*sigmab); normalizstion term kindof
%             bil_mask_amb = exp(double(-1*((amb_diffmask/sigmab1).^2)/(2*(sigmab1^2))));
            bil_mask_flash = exp(double(-1*((flash_diffmask).^2)/(2*(sigmab1^2))));      %/(sqrt(2*pi)*sigmab); normalizstion term kindof
            bil_mask_amb = exp(double(-1*((amb_diffmask).^2)/(2*(sigmab1^2))));
            filt_mask_flash = bil_mask_flash.*gauss_mask;
            norm_term_flash = sum(filt_mask_flash(:));
            filt_mask_amb = bil_mask_amb.*gauss_mask;
            norm_term_amb = sum(filt_mask_amb(:));
            Ajoint_mask = (double(amb_mask).*filt_mask_flash)/norm_term_flash;
            Abase_mask = (double(amb_mask).*filt_mask_amb)/norm_term_amb;
            Fbase_mask = (double(flash_mask).*filt_mask_flash)/norm_term_flash;
            Ajoint(i-bias,j-bias,:) = sum(sum(Ajoint_mask,1),2);
            Abase(i-bias,j-bias,:) = sum(sum(Abase_mask,1),2);
            Fbase(i-bias,j-bias,:) = sum(sum(Fbase_mask,1),2);
            if i ==60 && j ==60
                xxx=1;
            end
        end
        %i
    end
end

