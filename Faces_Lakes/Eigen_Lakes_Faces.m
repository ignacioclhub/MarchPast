

clear;
close all;

%// Read in image - make sure you convert it to double
%%%% Elija (1/0) paisaje (montania o lago)  o face (pelada o menton)
%%your choice
ver_face=0;
ver_pelada=0;
ver_montania=1;
%%//---------------------------------------------------------------

if ver_face == 0
B1 = imread('Lago-Moraine.jpeg'); 
B1 = rgb2gray(B1);
        if ver_montania == 1
    B=zscore(double(B1(1:200,:)')); %montania
        else
    B=zscore(double(B1(end-200:end,:)'));%lago
        end

else

B1 = imread('myface.tiff');
B1 = rgb2gray(B1);
            if ver_pelada == 1
            B=zscore(double(B1(1:200,:)')); %frente
            else
            B=zscore(double(B1(end-200:end,:)'));%menton
            end
 
end

[nn,mm]=size(B)

figure(1)
imshow(B')
xlabel('Fake Time axis')
ylabel('Fake ROIs')

  
 [T,N]=size(B);
 sd=std(B(:));
% Compute the boundaries from Marchenko Pastur 
c=N/T; 
a=(sd^2)*(1-sqrt(c))^2;
b=(sd^2)*(1+sqrt(c))^2;
%// Calculate covariance matrix
sigma = cov(B);
%// Compute the eigenvalues and eigenvectors of the covariance matrix
[A,D,W] = eig(sigma);
vals = diag(D);

figure
subplot(3,1,1)
pcolor(sigma)
xlabel('pixel I')
ylabel('pixel J')

shading flat
colorbar
title('Cov. Matrix')

%%// Compute the IPR
    for i=1:size(A,2)
     IPR(i)=0;
     for j=1:N
         IPR(i)=IPR(i)+A(j,i)^4;
     end
      
    end
   
     
 
 mP=find(vals > a & vals < b);  % encuentro los eigenvalues que son ruido (dentro del MP predictions)
subplot(3,1,3)
loglog(vals,IPR,'x',vals(mP),IPR(mP),'or')  % los que son rudio les pongo un circulo rojo
xlabel('\lambda')
ylabel('IPR')
title('IPR versus \lambda (red=Marchenko Pastur)')
%// Sort the eigenvalues
[sorteo,ind] = sort(abs(vals), 'descend');

%% Plot the rank de eigenvalues

subplot(3,1,2)
loglog(sorteo,'-ob')
xlabel('Rank')
ylabel('\lambda')
hold on
title('Ranked \lambda (red= Marchenko Pastur)')
%%  

 
 %% -----Tambien se puede computar lo mismo usando SVD, just for comparison (que debe dar exactamente lo mismo)
[U1,S1,V1] = svd(sigma);  
vals2=diag(S1);  %%
[sorteo2,ind2] = sort(abs(vals2), 'descend');
subplot(3,1,2)
loglog(sorteo2,'xg')  %% ploteo los eigenvalues computados usando el SVD
%%  those from Marchenko Pastur with red
mP=find(sorteo > a & sorteo < b);
loglog(mP,sorteo(mP),'<r')
hold on
grid on
 
%%// Recnstructing the time series (aka the pictire) using different number of eigenvalues/eigenvectors 

%// Rearrange eigenvectors usando el indice
Asort = A(:,ind);

%% If you wish to  see the effects of randomzing the eigenvectors 
%% over the reconstruction of the image uncoment the next line (shuffle) 
% Asort= shuffle(Asort,2);
%//  Mean subtracted data
Bm = bsxfun(@minus, B, mean(B,1));  
%// Reproject data onto principal components
Bproject = Bm*Asort;
%%%%%%%// Begin reconstruction 

counter = 1; 
figure;   
for k = [ 1 2 3 4 5 max(ind) ]  % Here choose wich k eigenvectors 
    %// Extract out highest k eigenvectors
    Aq = Asort(:,1:k);
    %// Project back onto the original domain
    picture_out = bsxfun(@plus, Bproject(:,1:k)*Aq.', mean(B, 1));
    %// Place projection onto right slot and show the image
    subplot(3, 2, counter);
    counter = counter + 1;
    imshow(picture_out',[]);  %%% ojo lo rotamos solo para graficar  y ver el paisaje o la face
    %imshow(picture_out')
    xlabel('fake time ->')
    ylabel('fake Rois ->')
    title(['k = 1 to ' num2str(k)]);
end

%% Proyectar los dos principal components, just for fun


Aq1 = Asort(:,1);  % o los que quieras 
Aq2 = Asort(:,2);
 
 figure
  plot(Aq1,Aq2,'-x')
 xlabel('Eigenvector 1')
  ylabel('eivenvector 2')
%  

%  for j=1:201
%      
%    plot(Asort(:,j)+(j/20),'-k');
% axis tight
% hold on
%  end
%  
