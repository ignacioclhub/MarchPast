clear
close all

%% config 

%Load subjects, change % to change groups.
    load AD10.mat; schaeferSignals=ts; %10 AD concatenados
    %load CN10.mat; schaeferSignals=ts; %10 CN concatenados
    %load MCI10.mat; schaeferSignals=ts; %10 MCI concatenados

% to identify subjects
load Schaefer2018_1000Parcels_17Networks_order_FSLMNI152_2mm.Centroid_RAS.mat;
Combined = strcat(ROI_Label, "_", string(Hemisphere));
RSN500=RSN_Label(1:500);
rsn_uniq=unique(RSN_Label); %nombre de las RSN para plotear

for rs=1:size(rsn_uniq,1) %saves matrix RSNxT just in case you want to work only with Resting State Networks
    schaRSN(rs,:)=mean(schaeferSignals(:,find(RSN_Label==char(rsn_uniq(rs))))');
end

%Script uses only 500 of 1000 ROIs, change it as desired
x=schaeferSignals(:,1:500); 
%x=schaRSN; %uncomment to work with RSN

%% Bucle con datos reales (figure 1) y shifteados (figure 2)
for prueba=1:2  % 1 is raw data, 2 is shufled data
    if prueba == 2  % aqui las shufleas
        for n=1:N
            shift=randi(T);
            x(n,:)=circshift(x(n,:),shift);
        end
        x=x';
    end

    [x, r, l, A] = preparemarchenko(x); %This function computes eigenvalues, vectors, IPR...
  
    s=std(x(:));
    [N,T]=size(x); %N es el numero de ROis y T el de samples
    c=N/T;


    % Probability Density Function
    n=50;     % number of points for measurement.
    a=(s^2)*(1-sqrt(c))^2;     % Boundaries -
    b=(s^2)*(1+sqrt(c))^2;    % Boundaries +
    [f,lambda]=hist(l,linspace(a,b,n));
    f=f/sum(f);     % Normalization la densidad

    % Theoretical pdf evaluated in the n points between lambda -a y lambda + b
    ft=@(lambda,a,b,c) (1./(2*pi*lambda*c*s^(2))).*sqrt((b-lambda).*(lambda-a));
    F=ft(lambda,a,b,c);
    % Processing numerical pdf
    F=F/sum(F); %Normalization
    F(isnan(F))=0;
    
%% Results
    figure(prueba)

    subplot(221)
        ll=sort(l,'descend');
        loglog(ll,'-.')
        hold on
        arriba=find(ll > b);
        if prueba == 1
            loglog(ll(arriba),'-og')  ;% derecha de MPasture
        else
        loglog(ll(arriba),'-or')
        end
        length(arriba)
        grid on
        hold on
        axis tight
    subplot(222)
        if prueba == 1
            h=bar(lambda,f);
            set(h,'FaceColor',[.75 .75 .8]);
            set(h,'LineWidth',0.25);
            xlabel('Eigenvalue \lambda');
            ylabel(' Probability Density Function f(\lambda)');
            hold on;
            plot(lambda,F,'g','LineWidth',2);
            hold on;
            title('Raw')
            lmin=min(l);
            lmax=max(l);
            axis([-1 1*lmax 0 max(f)+max(f)/4]);
        else
            h=bar(lambda,f);
            hold on;
            set(h,'FaceColor',[.5 .5 .8]);
            set(h,'LineWidth',0.5);
            plot(lambda,F,'r','LineWidth',2);
            lmin=min(l);
            lmax=max(l);
            axis([-1 1*lmax 0 max(f)+max(f)/4]);
            title('Shuffle')
        end
    subplot(223)
            %sortea IPR para buscar ROIs
            clear IPR L
            for i=1:size(A,2)
                IPR(i,1)=0;
                IPR_roi(i,1)=0;
                for j=1:N
                    IPR_roi(i,1)=IPR_roi(i,1)+A(i,j)^4; % de ROIs
                  IPR(i,1)=IPR(i,1)+A(j,i)^4; %original
                end
            end
            [IPRsort, index]=sort(IPR);
            IPRrois=IPR;
            IPRroi_label=ROI_Label(index);
            IPRrsn_label=RSN_Label(index);
            num_plots=17;
            ax=1:1:500;
            colors = jet(num_plots); % 'jet' es solo un ejemplo de mapa de colores. Puedes usar otros como 'hsv', 'hot', 'cool', etc.
        if prueba == 1
            for i=1:17
                plot(ax(find(RSN500==char(rsn_uniq(i)))),IPR_roi(find(RSN500==char(rsn_uniq(i)))),'o','Color', colors(i, :)); hold on
            end
            ylim([0,.2])
            set(gca,'YScale','log')
            ylabel('IPR ROI')
            legend(rsn_uniq)
        end
        if prueba == 2
            for i=1:17
                plot(ax(find(RSN500==char(rsn_uniq(i)))),IPR_roi(find(RSN500==char(rsn_uniq(i)))),'o','Color', colors(i, :)); hold on
            end
            set(gca,'YScale','log')
            ylim([0,.2])
            ylabel('IPR')
            legend(rsn_uniq)
        end
        subplot(224)
                %suma rangos.
                for i=1:17
                    RSNsumIPR(i)=sum(IPR_roi(RSN500==rsn_uniq(i)));
                    ranksum(i)=sum(index(IPRrsn_label==rsn_uniq(i)))/size(find(RSN500==rsn_uniq(i)),1);
                end
                [RSNsumIPR_sort, index] = sort(RSNsumIPR,'descend');
                plot(RSNsumIPR_sort,'o-','MarkerFaceColor','auto')                
             %   ylim([0,2])
                xticks([1:1:17]);
                xticklabels(rsn_uniq(index));
                ylabel('IPR sum');
                xlabel('RSN')



end