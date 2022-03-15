%% ICoRD v.1.1

clear all
close all
clc

% Initial parameters 
Saving_name="Final_result_name";
NN=10;              % The number of iterations. it depends on SNR of target neuron
Fs=10;                % Sampling frequency
tau_d = 0.550;     % decaying tau
tau_r = 0.179;       % rising tau
Early_stop=false;  
 
tic;
%% 1. Data load

% .mat or .tif can be loaded
[file,path] = uigetfile('*.mat;*.tif',...
               'Select an calcium imaging movie');

ismat=strfind(file,'.mat');
           if isempty(ismat)==1 %tif
                Data=imread([path,file]);% tif
                Y=Data;
                if size(Data,3)==1
                    initial_zero=5000; % Inital frame
                    numPlanes=2999; % The number of frames to analyze
                    for z = initial_zero:initial_zero+numPlanes
                        Data(:,:,z-initial_zero+1)=imread([path,file],z);
                    end
                      Y=double(Data);
                end
           else
                 Data=load([path,file]);%mat
                 Y=double(Data.Y);
           end
        
im2.sY=size(Y); 

      corr_image=zeros(im2.sY(1)-1,im2.sY(2)-1);
   for d=1+1:im2.sY(1)-1
     for f=1+1:im2.sY(2)-1
         q_1=corrcoef(reshape(Y(d,f,:),[1 im2.sY(3)]),reshape(Y(d-1,f,:),[1 im2.sY(3)]));
         q_2=corrcoef(reshape(Y(d,f,:),[1 im2.sY(3)]),reshape(Y(d+1,f,:),[1 im2.sY(3)]));
         q_3=corrcoef(reshape(Y(d,f,:),[1 im2.sY(3)]),reshape(Y(d,f+1,:),[1 im2.sY(3)]));
         q_4=corrcoef(reshape(Y(d,f,:),[1 im2.sY(3)]),reshape(Y(d-1,f-1,:),[1 im2.sY(3)]));
         corr_image(d,f)=mean([q_1(1,2) q_2(1,2) q_3(1,2) q_4(1,2)]);
     end
 end

im2.Cn=max(Y,[],3); 
mean_image=mean(Y,3);

figure(1);subplot(1,3,1);imagesc(im2.Cn);title('Max image');
subplot(1,3,2);imagesc(mean_image);title('Mean image');
subplot(1,3,3);imagesc(corr_image);title('Correlation image');
  
im2.thresholds=0:0.01:1-0.02;  
im2.SNR=zeros(size(im2.thresholds,2),1);
im2.optimal_threshold=0;
im2.NumData=0;

Y2d = reshape(Y,[im2.sY(1)*im2.sY(2) im2.sY(3)]); %transform to 2D
CN1d = reshape(im2.Cn,[im2.sY(1)*im2.sY(2) 1]); %transform to 1D

%% Correlation calculation and thresholding

n=0;hold on;button=1; 
%Mouse clicking
    while button==1,
      figure(2);gg=imagesc(mean_image);   gg.AlphaData=.95; hold on;
      title({'1. Select represntative pixel of target neuron with your mouse click (Right mouse: last neuron)'; '2. Square selection (left mouse drag & double click)'})
       [x,y,button]=ginput(1);               
       plot(x,y,'mo');

       c = strcat(int2str(x),' , ', int2str(y),' , ', int2str(n+1));
       text(x + 1,y,c); 

       %ROI pointer save
       n=n+1;  
       im2.Pointer(n,1)=n;
       im2.Pointer(n,2)=x;
       im2.Pointer(n,3)=y;

       % for ROI Detect
       h= imrect;  
       im2.ROIs(n,:)=wait(h);   %square save
       
         %% Draw pixels to reject
          figure(2); %imagesc(mean_image);
          title('3. Draw the region to reject (Drag to make circle)');
               h= imellipse;  
               Reject=wait(h);   %square save
        test_reject(:,:,n)=mean_image;
        for i=1:size(Reject,1)
                test_reject(round(Reject(i,2)),round(Reject(i,1)),n)=10;
        end
         reject_x=round(mean(Reject(:,2))); %y axis location
        b=round(max(Reject(:,2)))-round(min(Reject(:,2))); %y axis length
         half_b=b/2;
        reject_y= round(mean(Reject(:,1))); %x axis location
        a=round(max(Reject(:,1)))-round(min(Reject(:,1))); %x axis length
          half_a=a/2;
          clear Reject
         %reject pixels
        test_reject(:,:,n)=mean_image;
        for x=1:size(test_reject,1)
            for y=1:size(test_reject,2)
                if ((y-reject_y)^2)/half_a^2+((x-reject_x)^2)/half_b^2<=1
                test_reject(x,y,n)=10;
                end
            end
        end
    end
%         figure();imagesc(test_reject(:,:,1))
%     figure();imagesc(test_reject(:,:,2))

 for nn=1:size(im2.ROIs,1)
ROI_x1(nn) =im2.ROIs(nn,1);
ROI_y1(nn) =im2.ROIs(nn,2);
ROI_x2(nn)=ROI_x1(nn)+im2.ROIs(nn,3);
ROI_y2(nn)=ROI_y1(nn)+im2.ROIs(nn,4);
 end
 NumofROIs=size(im2.Pointer,1);


    im2.ROIY_signal = zeros(NumofROIs,im2.sY(3));
    im2.ROIY_noise = zeros(NumofROIs,im2.sY(3));
    im2.ROIYY = zeros(NumofROIs,im2.sY(3));
    im2.ROIYY_signal = zeros(NumofROIs,im2.sY(3));
    im2.ROIYY_noise = zeros(NumofROIs,im2.sY(3));
    im2.ROINN = zeros(NumofROIs,im2.sY(3));
    im2.ROIY = zeros(NumofROIs,im2.sY(3));
    im2.ROIN = zeros(NumofROIs,im2.sY(3));
    im2.ROINN = zeros(NumofROIs,im2.sY(3));
    

for ROI_i=1:NumofROIs

im2.Cn1=im2.Cn;    
refD_temp = reshape(Y(round(im2.Pointer(ROI_i,3)),round(im2.Pointer(ROI_i,2)),:),[1,im2.sY(3)]);
refD_temp_mean=mean(mean(Y(round(im2.Pointer(ROI_i,3))-1:round(im2.Pointer(ROI_i,3))+1,round(im2.Pointer(ROI_i,2))-1:round(im2.Pointer(ROI_i,2))+1,:),1),2);
refD_temp_mean=reshape(refD_temp_mean,[1,im2.sY(3)]);
refD_temp=refD_temp_mean;%average of around 9 pixels to calculate first reference

%% tau
dt = 1/Fs;
A = [-(2/tau_d+1/tau_r), - (tau_r+tau_d)/(tau_r*tau_d^2); 1 0];
lc = eig(A*dt);
ld = exp(lc);
g2 = [sum(ld),-prod(ld)];
h2 = (1-exp(-dt/tau_r))*exp(-dt/tau_d);

[refD,cb,c1,~,~,im.Espikes] = constrained_foopsi(refD_temp',[],[],g2);    
Iter_results=zeros(im2.sY(3),NN);
Iter_results_raw=zeros(im2.sY(3),NN);
Iter_results_CD=zeros(im2.sY(3),NN);

time_1=0:1/Fs:size(refD,1)/Fs -1/Fs;
figure();plot (time_1,normalize(refD_temp));hold on;plot(time_1,normalize(refD),'LineWidth',1.5); 
legend('Raw of initial reference', 'CD for initial reference ')
xlabel('Time(s)','FontSize',12); ylabel('Normalized amplitude','FontSize',12)

%ROI detection according to the thresholds
th_ROI=zeros(size(im2.thresholds,2),im2.sY(1),im2.sY(2));    
im2.Detected_ROI=zeros(size(im2.thresholds,2),im2.sY(1),im2.sY(2));

 for th = 1:size(im2.thresholds,2)
     im2.Cn1=im2.Cn;
       
       h=1;xx=0;yy=0;p=0;q=0; 
        for yy = round(im2.ROIs(1,2)):round(im2.ROIs(1,2)+im2.ROIs(1,4))% y axis
            for xx =round(im2.ROIs(1,1)):round(im2.ROIs(1,1)+im2.ROIs(1,3))  % x axis
                if test_reject(yy,xx,ROI_i)==10
                    continue
                end
                 comD = reshape(Y(yy,xx,:),[1,im2.sY(3)]);
                 corND_1 = corrcoef(refD,comD); % refD_temp can be used
                 Corr(1,h) = corND_1(1,2);

                  if Corr(1,h)>im2.thresholds(th) %thresholing
                     p=p+1;
                     im2.Cn1(yy,xx)=20000;%1.3
                     im2.ROIY(1,:)=im2.ROIY(1,:)+comD;
                  else 
                     q=q+1;
                     im2.ROIN(1,:)=im2.ROIN(1,:)+comD;
                  end
                  h=h+1;
           end
        end
        th_ROI(th,:,:)=im2.Cn1;
        im2.NumData(1,th)=p;
        
    im2.ROIYY(th,:) = im2.ROIY(1,:)/p;
    im2.ROINN(th,:) = im2.ROIN(1,:)/q;
    im2.ROIY = zeros(NumofROIs,im2.sY(3));
    im2.ROIN = zeros(NumofROIs,im2.sY(3));
 end

assignin('base','im2',im2);
disp('ROI selection is done');

im.numFrames=1;
im.numPlanes=im2.sY(3);


tau_rise1 = 0.179;
tau_decay1= 0.550;
dt = 1/Fs;
[g2,h2] = tau_c2d(tau_rise1,tau_decay1,dt);
Ideal=zeros(im2.sY(3),size(im2.thresholds,2));
References=zeros(size(im2.thresholds,2),im2.sY(3));
CD_each_N=zeros(size(im2.thresholds,2),im2.sY(3));
Average_each_N=zeros(size(im2.thresholds,2),im2.sY(3));
First_raw=zeros(size(im2.thresholds,2),im2.sY(3));


 for th=1:size(im2.thresholds,2)
 im.NdF = [];
%
if isnan(im2.ROIYY(th,1))
break;
end

%normalization
im2.ROID = im2.ROIYY(th,:);
im.NdF = im2.ROID-(min(im2.ROID)); 
im.NdF = im.NdF/max(im.NdF);

%foopsi algorithm
[im.filtPix,cb,c1,~,~,im.Espikes] = constrained_foopsi(im.NdF,[],[],g2);
 Raw_s=im.NdF';
Ideal_s=im.filtPix;
         if isnan(im.filtPix(1,1))==1
             continue;
         end
% im2.filtPix(th,:)=im.filtPix;
        Ideal(:,th)=im.filtPix;
        corr_temp=corrcoef(im.filtPix,im2.ROID ); % corr estimation
        im2.SNR(th)=corr_temp(1,2);

 im2.SNR=nonzeros(im2.SNR);
    if im2.NumData(1,th) ==1
        im2.SNR(th:size(im2.SNR,1))=[];
       
        break;
    end

Average_each_N(th,:)=Raw_s;
References(th,:)=Ideal_s;
Ideal_s=[];


 end
[c,d]=sort(nonzeros(im2.SNR),'descend');
if c(1)==0
    a=c(2);
    b=d(2);
else
    a=c(1);
    b=d(1);
end
 im2.optimal_threshold=im2.thresholds(b);


%% Iterations

im2.optimal_SNR=zeros(1,NN);
im2.optimal_th=zeros(1,NN);
im2.optimal_SNR=a;
im2.optimal_th=b;
% Iter_results=zeros(im2.sY(3),NN);
Iter_results_raw(:,1)=Average_each_N(im2.optimal_th,:);
Iter_results_CD(:,1)=References(im2.optimal_th,:);
Iter_results(:,1)=References(im2.optimal_th,:);

Iter_results_SNR=zeros(NN,size(im2.thresholds,2));
Iter_results_SNR=Iter_results_SNR;
Iter_results_SNR(1,1:size(im2.SNR,1))=im2.SNR;
Iter_results_ROI=zeros(NN,im2.sY(1),im2.sY(2)); 
Iter_results_ROI(1,:,:)=th_ROI(im2.optimal_th,:,:);
Iter_results_NumData=zeros(NN,size(im2.NumData,2));
Iter_results_NumData(1,:)=im2.NumData;
Iter_spikes=zeros(im2.sY(3),NN);
Iter_spikes(:,1)=im.Espikes;

information_difference=zeros(2,NN);
P_val=zeros(2,NN);
         window_size=5;
  mean_temp=zeros(NN-window_size,1);
for N=2:NN
    N
new_ref=Iter_results(:,N-1);%new reference
% deconvolution_sig=[];
im2.Detected_ROI=zeros(size(im2.thresholds,2),im2.sY(1),im2.sY(2));
 for th = 1:size(im2.thresholds,2)
     im2.Cn1=im2.Cn;
       refD = new_ref';
       h=1;xx=0;yy=0;p=0;q=0; %h :Column of the picked coordinates (Data)
        for yy = round(im2.ROIs(1,2)):round(im2.ROIs(1,2)+im2.ROIs(1,4))%y-coordinate shift variable
            for xx =round(im2.ROIs(1,1)):round(im2.ROIs(1,1)+im2.ROIs(1,3))  % x-coordinate shift variable
                if test_reject(yy,xx,ROI_i)==10
                    continue
                end
                comD = reshape(Y(yy,xx,:),[1,im2.sY(3)]);
                 corND_1 = corrcoef(refD,comD);
                 Corr(1,h) = corND_1(1,2);

                  if Corr(1,h)>im2.thresholds(th) %thresholing
                     p=p+1;
                     
                     im2.Cn1(yy,xx)=20000;
                     im2.ROIY(1,:)=im2.ROIY(1,:)+comD;
                  else 
                     q=q+1;
                     im2.ROIN(1,:)=im2.ROIN(1,:)+comD;
                  end
                  h=h+1;
           end
        end
        th_ROI(th,:,:)=im2.Cn1;
        im2.NumData(1,th)=p;
        
        
    im2.ROIYY(th,:) = im2.ROIY(1,:)/p;
    im2.ROINN(th,:) = im2.ROIN(1,:)/q;
    im2.ROIY = zeros(NumofROIs,im2.sY(3));
    im2.ROIN = zeros(NumofROIs,im2.sY(3));

 end


         for th=1:size(im2.thresholds,2)
         im.NdF = [];
        %
        %baseline subtraction
        im2.ROID = im2.ROIYY(th,:);
         if isnan(im2.ROID(1,1))==1;
             break
         end

        %normalization
        im.NdF = normalize(im2.ROID);

        %foopsi algorithm
        % th
         Raw_s=im.NdF;
        Ideal_s=new_ref;
        % im2.filtPix(th,:)=im.filtPix;



        if isnan(im2.ROID(1,1))~=1;
        [im.filtPix,cb,c1,~,~,im.Espikes] = constrained_foopsi(Raw_s,[],[],g2);
        end
         
         if isnan(im.filtPix(1,1))==1
             continue;
         end

        corr_temp=corrcoef(im.filtPix,im2.ROID ); 
        im2.SNR(th)=corr_temp(1,2);
          
        CD_each_N(th,:)=normalize(im.filtPix');
        Average_each_N(th,:)=im.NdF ; 
            
    im2.SNR=nonzeros(im2.SNR);
    if im2.NumData(1,th) ==1
        im2.SNR(th:size(im2.SNR,1))=[]; 
        break;
    end
        Ideal_s=[];

         end

Iter_results_SNR(N,1:size(im2.SNR,1))=im2.SNR;
Iter_results_SNR(N,1:size(im2.SNR,1))=nonzeros(Iter_results_SNR(N,:))';
[c,d]=sort(nonzeros(im2.SNR),'descend');
if c(1)==0
    a=c(2);
    b=d(2);
else
    a=c(1);
    b=d(1);
end

for ggg=1:size(size(im2.SNR,1))
if Iter_results_SNR(N,ggg)==0
    temp(ggg)=Iter_results_SNR(N,ggg);
    temp(ggg)=[];    
    Iter_results_SNR(N,ggg)= temp(ggg);   
end
end

im2.optimal_threshold=im2.thresholds(b);
im2.optimal_SNR(:,N)=a;
im2.optimal_th(:,N)=b;


Iter_results_raw(:,N)=Average_each_N(im2.optimal_th(:,N),:);
Iter_results_CD(:,N)=CD_each_N(im2.optimal_th(:,N),:);


corr_temp=corrcoef(Iter_results_raw(:,N),Iter_results_raw(:,N-1));
 P_val(1,N)=corr_temp(1,2);
information_difference(1,N)=1-corr_temp(1,2);



 P_val(2,N)=erfc(abs(atanh( P_val(1,N))-atanh( P_val(1,N-1)))/(sqrt(2)*sqrt(2/497)));

threshold_for_update_decision=0.02;
    if information_difference(1,N)<threshold_for_update_decision|| information_difference(1,N)>information_difference(1,N-1) %When the ID is saturated,  change the method!!

    if Iter_results(:,N-1)==Iter_results_raw(:,N-1)
    Iter_results(:,N)=Iter_results_CD(:,N); 
    information_difference(2,N)=1;
    else
    Iter_results(:,N)=Iter_results_raw(:,N);       
        information_difference(2,N)=0;
    end
else %When the ID is  still not saturated,  keep the method!!
    if Iter_results(:,N-1)==Iter_results_raw(:,N-1)
    Iter_results(:,N)=Iter_results_raw(:,N);       
            information_difference(2,N)=0;
    else
    Iter_results(:,N)=Iter_results_CD(:,N); 
        information_difference(2,N)=1;
    end
end

Iter_results_ROI(N,:,:)=th_ROI(im2.optimal_th(:,N),:,:);
Iter_results_NumData(N,:)=im2.NumData(im2.optimal_th(:,N));
Iter_spikes(:,N)=im.Espikes;
      if (im2.optimal_SNR(:,N)==Inf)
            break
      end

      
      if N>window_size-1 %From 4, we update the estimated corr of the moving window and decide when to stop.
          
         reassigned_estimated_SNR_temp=zeros(1,window_size);
         temp_min=min(Iter_results_raw(:,N));
         temp_max=max(Iter_results_raw(:,N)-temp_min);
         [x2,cb,c1,~,~,im.Espikes] = constrained_foopsi(normalize(Iter_results_raw(:,N)),[],[],g2);
         iii=0;
            for d=N-(window_size-1):N
                  iii=iii+1;  
         corr_temp=corrcoef(x2,Iter_results_CD(:,N) ); % corr estimation
         reassigned_estimated_SNR_temp(iii)=corr_temp(1,2);
            end
                 mean_temp(N)=mean(reassigned_estimated_SNR_temp);

%          %% Criteria for stopping iteration
%          if Early_stop==true
%                      if N>window_size 
%                            if mean_temp(N-1)>mean_temp(N)
%                                NN=N;
%                                break;
%                            end
%                      end
%          end
    end
end

figure();plot(Iter_results_NumData(:,1),'-o');title('Num of ROIs in iteration');xlabel('The number of iterations','FontSize',12);ylabel('The number of pixels','FontSize',12);

%ROI trends
figure();
for(i=1:NN)
   num_temp= floor(NN/5);
   if num_temp==0
       num_temp=1;
   end
    subplot(num_temp+1,5,i);imagesc(reshape(Iter_results_ROI(i,:,:),[im2.sY(1) im2.sY(2)])); title(['Iteration ',num2str(i)]);
%     xlim([ROI_x1 ROI_x2])
%     ylim([ROI_y1 ROI_y2])
end

%% Calculate and update im2.SNR with a new reference
 reassigned_estimated_SNR=zeros(1,NN);
 temp_min=min(Iter_results_raw(:,NN));
 temp_max=max(Iter_results_raw(:,NN)-temp_min);
 [x2,cb,c1,~,~,im.Espikes] = constrained_foopsi(normalize(Iter_results_raw(:,NN)),[],[],g2);
for N=1:NN      
        corr_temp=corrcoef(x2 ,Iter_results_CD(:,N)); % SNR estimation 
        reassigned_estimated_SNR(N)=corr_temp(1,2);    
end

[a,b]=max(reassigned_estimated_SNR);
% reassigned_estimated_SNR
figure();plot(im2.optimal_SNR,'-o'); hold on;
plot(reassigned_estimated_SNR,'-o');xlabel('The number of iterations','FontSize',12); ylabel('SNR (dB)','FontSize',12)
legend('Estimation corr coef','reassigned corr coef');


%   figure();
  std_temp=zeros(NN-window_size,1);
 for k=window_size:NN
   std_temp(k)= mean(reassigned_estimated_SNR(k-window_size+1:k)); 
 end
 std_temp(1:window_size-1)=NaN;
% plot(std_temp);




%% Figures
  [Max_SNR opti_Iter ]=max(reassigned_estimated_SNR)


figure();plot(time_1,normalize(refD_temp),'Linewidth',1.5,'Color',[0.7 0.7 0.7]); hold on;
plot(time_1,normalize(Iter_results_CD(:,1)),'Linewidth',1.5,'Color',[0 0.54 0.77]);
hold on; plot(time_1,normalize(Iter_results_CD(:,opti_Iter)),'Linewidth',1.5,'Color',[1 0.6 0]);legend('First raw','First deconvolution', 'Final result')
xlabel('Time (s)','FontSize',12); ylabel('Normalized amplitude','FontSize',12);title("Neuron"+num2str(ROI_i))


figure();plot(time_1,normalize(Iter_results_raw(:,opti_Iter)),'Linewidth',1.5,'Color',[0 0.54 0.77]);hold on; 
plot(time_1,normalize(Iter_results_CD(:,opti_Iter)),'Linewidth',1.5,'Color',[1 0.6 0]);legend('Final raw', 'Final result')
 xlabel('Time (s)','FontSize',12); ylabel('Normalized amplitude','FontSize',12);title("Neuron"+num2str(ROI_i))

assignin('base','im2',im2);
assignin('base','im',im);



 figure();
 for i=1:NN

 plot(im2.thresholds(1:size(nonzeros(Iter_results_SNR(i,:)))),nonzeros(Iter_results_SNR(i,:))); hold on;
%  [aa,bb]=max(Iter_results_SNR(i,:));
aa=im2.optimal_SNR(:,i);
bb=im2.optimal_th(:,i);

 plot(im2.thresholds(bb),aa,'r*');
c = strcat('(',num2str(im2.thresholds(bb)),' , ', num2str(aa),')_*',num2str(i));
text(im2.thresholds(bb)+0.05,aa,c); 
end
 xlabel('Threshods','FontSize',12); ylabel('SNR (dB)','FontSize',12)
 title('SNR(estimation) according to correlation thresholds in iterations');

toc;

Final.opti_Iter(ROI_i)=opti_Iter-1;
Final.final_proposed_temporal_raw(:,ROI_i)=normalize(Iter_results_raw(:,opti_Iter-1));
Final.final_proposed_temporal_CD(:,ROI_i)=normalize(Iter_results_CD(:,opti_Iter-1));
Final.iter_ROIs(:,:,:,ROI_i)=Iter_results_ROI;
Final.Background(ROI_i,:)=normalize(refD_temp);
Final.final_proposed_spikes(:,ROI_i)=normalize(Iter_spikes(:,NN));
Final.final_estimatedSNR(ROI_i,:)=reassigned_estimated_SNR;
Final.Iter_results_SNR(:,:,ROI_i)=Iter_results_SNR;
Final.information_difference(:,:,ROI_i)=information_difference;
Final.threshold_for_update_decision(ROI_i)=threshold_for_update_decision;
Final.P_val(ROI_i)=P_val(1);
end

Final.Max_image=max(Y,[],3);
Final.Mean_image=mean_image;
Final.final_time=time_1;
Final.regionofcalculation=[ROI_x1,ROI_y1,ROI_x2,ROI_y2];
Final.iter_num=1:NN;
Final.thresholds=im2.thresholds;
save(char(Saving_name),'Final');

 for ROI_i=1:NumofROIs
figure(); title("Neuron"+num2str(ROI_i))
for(i=1:NN)
       num_temp= floor(NN/5);
   if num_temp==0
       num_temp=1;
   end
    subplot(num_temp,5,i);imagesc(reshape(Final.iter_ROIs(i,:,:,ROI_i),[im2.sY(1) im2.sY(2)])); title(['Iteration ',num2str(i)]);
%     xlim([ROI_x1(ROI_i) ROI_x2(ROI_i)])
%     ylim([ROI_y1(ROI_i) ROI_y2(ROI_i)])
end
 end
 
 for ROI_i=1:NumofROIs
 figure(); 
imagesc(reshape(Final.iter_ROIs(i,:,:,ROI_i),[im2.sY(1) im2.sY(2)])); title("ROI of neuron"+num2str(ROI_i));
end


