
% This function gets fmri timeseries and transform them to get eigenvalues
% , eigenvectors and IPR(inverse participation index)
% INPUT: Matrix t (TR of timeseries) x R (ROIS)


function [signals, rmatrix, eigenvalues, eigenvectors, IPR] = preparemarchenko(x)

[T,N]=size(x); %tamaño de la matriz a trabajar
x=x'; %traspone la señal
for n=1:N %sustrae la media a las señales
        x(n,:)=x(n,:)-mean(x(n,:));
end
x=zscore(x); %normaliza la señal

[N,T]=size(x); %N es el numero de ROis y T el de samples
c=N/T; % Ratio of matrix dimensions
s=std(x(:)); %desviación de la señal
%trabajo con la matriz de covarianza
r=x*x'/T; % spectral matrix, covarianza de la señal.
r=(r+r')/2;
%r(1:1+size(r,1):end) = 0;   % put zero en la diagonal
r(isnan(r)) = 0;
r(isinf(r)) = 0;

[A,D] = eig(r); %compute eigenvalues and eigenvectors
l = abs(diag(D));

%inverse participation index
for i=1:size(A,2)
    IPR(i,1)=0;
    for j=1:N
      IPR(i,1)=IPR(i,1)+A(j,i)^4; 
    end
end





signals = x;
rmatrix = r;
eigenvalues = l;
eigenvectors = A;
IPR = IPR;