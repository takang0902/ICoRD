
%% Simulated neuron in calcium imaging

clear all;
close all
clc

SNR_in = -30; % Desired Signal-to-noise power ratio (dB) of simulated neuron 
Fs=10; % Sampling rate

trial = 10;
g_sigma = 10;
T=50*trial; %50 seconds
siz = [150, 150, T]; 
sigma_seed = 1;
num_neuron = 1;

tau_dd = 0.550; % tau d of GCaMP6 
tau_rr = 0.179; % tau r of GCaMP6 
zero_array = zeros(1,50);

K = siz(1)*siz(2);
numFrames = 50;
numPlanes = T;
stimFrame = 5;
condition = 1;
dt = 1/Fs;
TT=50;

DotS = '.';
colors=[1 1 1;1 0 0;0.5 0.2 0.3;0 0 0;0.1 1 0.1];
temporal_model = zeros(130,50);

range = 1;  

for i = 1 : range

tau_d = 0.550;
tau_r = 0.179;
%% tau
A = [-(2/tau_d+1/tau_r), - (tau_r+tau_d)/(tau_r*tau_d^2); 1 0];
lc = eig(A*dt);
ld = exp(lc);
g2 = [sum(ld),-prod(ld)];
h2 = (1-exp(-dt/tau_r))*exp(-dt/tau_d);

%% 
if stimFrame == 0
    msig = filter(h2,[1,-g2],[zeros(1, stimFrame-1) 1 zeros(1, numFrames-stimFrame-1)]);
else
    msig = filter(h2,[1,-g2],[zeros(1, stimFrame-1) 1 zeros(1, numFrames-stimFrame)]);
end
temp_msig = filter(h2,[1,-g2],[ 1 zeros(1, TT-1)]);
Nmsig = msig/(rms(msig)*numFrames);
temp_Nmsig = temp_msig/(rms(temp_msig)*TT);
sGLM = transpose([ones(1,numFrames)/numFrames; Nmsig]);

temporal_model(i,:) = Nmsig;
end

%% Calcium decaying model
AA = [-(2/tau_dd+1/tau_rr), - (tau_rr+tau_dd)/(tau_rr*tau_dd^2); 1 0];
lcc = eig(AA*dt);
ldd = exp(lcc);
g22 = [sum(ldd),-prod(ldd)];
h22 = (1-exp(-dt/tau_rr))*exp(-dt/tau_dd);

%% 
if stimFrame == 0
    msigg = filter(h22,[1,-g22],[zeros(1, stimFrame-1) 1 zeros(1, numFrames-stimFrame-1)]);
else
    msigg = filter(h22,[1,-g22],[zeros(1, stimFrame-1) 1 zeros(1, numFrames-stimFrame)]);
end
temp_msigg = filter(h22,[1,-g22],[ 1 zeros(1, TT-1)]);
Nmsigg = msigg/(rms(msigg)*numFrames);
sGLMm = transpose([ones(1,numFrames)/numFrames; Nmsigg]);
sGLM_invv = pinv(sGLMm);
%% spatial
siz_y = siz(1); %150
siz_x = siz(2); %150

location_of_neurons_x = zeros(num_neuron,1);
location_of_neurons_y = zeros(num_neuron,1);
for i=1:num_neuron
    location_of_neurons_x(i) = 75;
    location_of_neurons_y(i) = 75; %center of neuron
end

sigma_x = (1+ 0.1*sigma_seed)*g_sigma;
sigma_y = (1+ 0.1*sigma_seed)*g_sigma; %cell size

[x,y] = meshgrid(1:siz_x, 1:siz_y);
spatial = zeros(siz_y, siz_x, num_neuron);
for m=1:num_neuron
    init = ((x-location_of_neurons_x(m))/sigma_x(1)).^2 + ((y-location_of_neurons_y(m))/sigma_y(1)).^2;
    temp = (exp(-init/2))/(sigma_x(1)*sigma_y(1)*2*pi);
    temp(temp<0.00000016232463096)=0; % 3sigma
    spatial(:,:,m) = temp;
end

spat_temp = reshape(spatial, K,1);
correct = find(spat_temp);
temtemtem = size(correct);

%% temporal

rng('shuffle'); 
C_model = zeros(range, num_neuron, TT*trial);

for l = 1 : range

S = zeros(num_neuron, T); 
C = zeros(num_neuron, T); 
kernel = temporal_model(l,:);

spks_activity=zeros(1,500);  
        for n = 1:trial
           
            spks = zero_array;
            if n==1
            spks(5) = 0.4;
            spks(1) = 0.5;
            elseif n==3
            spks(1) = 0.7;
            spks(10) = 0.3; 
            elseif n==4
            spks(7) = 0.5;
            spks(20) = 0.3;
            spks(3) = 0.5; 
            elseif n==5
            spks(5) = 0.4;
            spks(1) = 0.7; 
            elseif n==9
            spks(1) = 0.4; 
            elseif n==7
            spks(9) = 0.4;
            spks(4) = 0.8; 
            elseif n==6
            spks(4)=1;
            else
            spks(1) = 0.56; 
            spks(15) = 0.3;
            end
            
           spks_activity(1+(n-1)*numFrames:n*numFrames) = spks;
            c = conv(spks,kernel);

            C(1,TT*(n-1)+1:TT*n) = c(1,1:TT);
        end 
C_model(l,:,:) = C;
end


%% spatio_temporal

spatial1=reshape(spatial,siz(1)*siz(2),[]);

simulation_result = zeros(range, siz(1), siz(2), T);
for o = 1:range
simulation_result_temp = spatial1*reshape(C_model(o,:,:),num_neuron,[]);
simulation_result(o,:,:,:) = reshape(simulation_result_temp, siz(1), siz(2), T);
end


simulation_result1 = reshape(simulation_result(1,:, :,:),siz(1),siz(2),siz(3));


%% power calculation
pow_result = simulation_result_temp;
pow_result = reshape(pow_result,siz(1),siz(2),500);
pow_temp = pow_result(1:150,1:150,:);
pow_temp_size = numel(pow_temp);

pow_sq = pow_temp.^2;
pow_sq_sum = sum(sum(sum(pow_sq)));
pow_pow = pow_sq_sum/pow_temp_size; %Mean power
pow_pow_sqrt = sqrt(pow_pow);
pow_db = mag2db(pow_pow_sqrt);
pow_noise = pow_db - SNR_in; %(dB) Mean noise power

%% add white Gaussian noise
noise_s = wgn(K,siz(3),pow_noise);
noise_3d = reshape(noise_s,siz(1),siz(2),siz(3));

add_noise = zeros(range, siz(1), siz(2), T);
for o = 1:range
for i = 1:T
    add_noise(o,:,:,i) = reshape(simulation_result(o,:,:,i),siz(1),siz(2),1) + noise_3d(:,:,i);
end
end
add_noise = add_noise + abs(min(min(min(min(add_noise)))));

noised_result = reshape(add_noise(1,:, :,:),siz(1),siz(2),siz(3));
time_1=0:1/Fs:size(noised_result,3)/Fs -1/Fs;

max_value=max(reshape(noised_result(75,75,:),1,siz(3)));
figure(); 
 plot(time_1,reshape(simulation_result1(75,75,:)/max_value,1,siz(3)),'Linewidth',0.8 ); % signal data
hold on; plot(time_1,reshape(noised_result(75,75,:)/max_value,1,siz(3)),'Linewidth',0.8); % raw data
hold on; plot(time_1,reshape(noised_result(75,75,:)/max_value,1,siz(3))-reshape(simulation_result1(75,75,:)/max_value,1,siz(3)),'Linewidth',0.8 ); % Noise data
legend('Singal','Signal + noise','Noise');
ylabel('Normalized amplitude','FontSize',12)
xlabel('Time(sec)','FontSize',12)

% Storage
savedata= noised_result; 
t=Tiff(char("_"+num2str(abs(SNR_in))+"dB.tif"),'w'); 
tagstruct.ImageLength = size(savedata,1);
tagstruct.ImageWidth = size(savedata,2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 64; %double
tagstruct.SampleFormat = 3; %double 
tagstruct.SamplesPerPixel = size(savedata,3);
tagstruct.RowsPerStrip = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software = 'MATLAB';
tagstruct % display tagstruct
setTag(t,tagstruct)
write(t,savedata);
close(t);

read_data=imread(char("_"+num2str(abs(SNR_in))+"dB.tif")); 

model_fig1=reshape(simulation_result1(75,75,:),1,siz(3));
      corr_image=zeros(size(noised_result,1)-1,size(noised_result,2)-1);
   for d=1+1:size(noised_result,1)-1
     for f=1+1:size(noised_result,2)-1
         q_1=corrcoef(reshape(noised_result(d,f,:),[1 size(noised_result,3)]),reshape(noised_result(d-1,f,:),[1 size(noised_result,3)]));
         q_2=corrcoef(reshape(noised_result(d,f,:),[1 size(noised_result,3)]),reshape(noised_result(d+1,f,:),[1 size(noised_result,3)]));
         q_3=corrcoef(reshape(noised_result(d,f,:),[1 size(noised_result,3)]),reshape(noised_result(d,f+1,:),[1 size(noised_result,3)]));
         q_4=corrcoef(reshape(noised_result(d,f,:),[1 size(noised_result,3)]),reshape(noised_result(d-1,f-1,:),[1 size(noised_result,3)]));
         corr_image(d,f)=mean([q_1(1,2) q_2(1,2) q_3(1,2) q_4(1,2)]);
     end
   end
   
figure();
subplot(1,3,1);imagesc(max(read_data,[],3));title('Max image');
subplot(1,3,2);imagesc(mean(read_data,3));title('Mean image');
subplot(1,3,3);imagesc(corr_image);title('Mean correlation image');
  

savedata2= savedata; 
savedata2=uint16(savedata2);%uint16
savedata2=savedata2(:,:,200);
t2=Tiff(char(num2str(SNR_in)+"dB_image.tif"),'w');
tagstruct.ImageLength = size(savedata2,1);
tagstruct.ImageWidth = size(savedata2,2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 16; %uint16
tagstruct.SampleFormat = Tiff.SampleFormat.UInt %uint16
tagstruct.SamplesPerPixel = size(savedata2,3);
tagstruct.RowsPerStrip = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software = 'MATLAB';
tagstruct % display tagstruct
setTag(t2,tagstruct)
write(t2,savedata2);
close(t2);

%% matfile storage
Y=savedata;
save("_"+char(num2str(abs(SNR_in))+"dB_image.mat"),'Y','time_1','model_fig1');





