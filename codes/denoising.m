clear all; close all; clc;
imflash = im2double(imread('C:\Users\Yang Shu\Desktop\520project\img\image1_flash.bmp'));
imambient = im2double(imread('C:\Users\Yang Shu\Desktop\520project\img\image1_noflash.bmp'));


Fr = imflash(:,:,1);
Fg = imflash(:,:,2);
Fb = imflash(:,:,3);
ambient_gradient = gradient(imambient);

[Ajointr, Abaser, Fbaser] = bilateral(Fr, imambient(:,:,1));%-ambient_gradient(:,:,1)*2);
[Ajointg, Abaseg, Fbaseg] = bilateral(Fg, imambient(:,:,2));%-ambient_gradient(:,:,2)*2);
[Ajointb, Abaseb, Fbaseb] = bilateral(Fb, imambient(:,:,3));%-ambient_gradient(:,:,3)*2);

Ajoint = cat(3,Ajointr,Ajointg,Ajointb);
Fbase = cat(3,Fbaser,Fbaseg,Fbaseb);
Abase = cat(3,Abaser,Abaseg,Abaseb);

eps = 0.02;
Fdetailr = (Fr+eps)./(Fbaser+eps);
Fdetailg = (Fg+eps)./(Fbaseg+eps);
Fdetailb = (Fb+eps)./(Fbaseb+eps);
Fdetail = cat(3,Fdetailr,Fdetailg,Fdetailb);
% Fdetail = Fdetail + gradient(imflash)*2;

shadowMask = shadowRem(imflash, imambient);
Ffin = (cat(3,(1-shadowMask),(1-shadowMask),(1-shadowMask)).*Ajoint.*Fdetail) + (cat(3,shadowMask,shadowMask,shadowMask).*Abase);
%Ffin = whiteBal(Ffin,imflash,2);

figure, imshow(imambient);
figure, imshow(imflash);
figure, imshow(Ffin);
