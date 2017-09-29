
xx=load('/home/pablo/local.git/pablo/res/Resultados_SART/Resultados_SART_4.txt');
tam=size(xx);
[minimum, indexmin]=min(xx(:));
[maximum2, indexmin2]=max(xx(:,3));
fprintf('El valor minimo es %d\n', minimum);
for i=1:tam(1)
  img2(int16(xx(i,1))+1, int16(xx(i,2))+1)=xx(i,3);  
    
   % if(xx(i,3)==0.00)
        
   %     xx(i,3)=(maximum2+minimum)/2; 
    %end
img(int16(xx(i,1))+1, int16(xx(i,2))+1)=xx(i,3);

    
end
[minimum, indexmin]=min(img2(:));
[maximum, indexmax]=max(img2(:));

tam2=size(img,1);

theta = 0:179;
N_theta = length(theta);
[R,xp] = radon(img,theta);

N1 = length(xp);
freqs=linspace(-1, 1, N1).';
myFilter = abs( freqs );
myFilter = repmat(myFilter, [1 N_theta]);

ft_R = fftshift(fft(R,[],1),1);
filteredProj = ft_R .* myFilter;
filteredProj = ifftshift(filteredProj,1);
ift_R = real(ifft(filteredProj,[],1));
filter_type='Ram-Lak';
I2 = iradon(R, theta, 'linear',filter_type, 1.0, tam2);
[minpad, indexpad]=min(I2(:));
[maxpad, indexpadmax]=max(I2(:));
title1=sprintf('Con Pad y filtro %s', filter_type);
imshow(I2, [minpad, maxpad]);title(title1);

figure;
imshow(img2,[minimum, maximum]);title('Sin Pad ni filtrado');
