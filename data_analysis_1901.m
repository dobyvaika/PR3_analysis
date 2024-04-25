function opto_data_analysis
% fail = -100;
mkdir results
output_folder = 'results\';
h = 0.02; Fs = 1/h;

    listing=dir;
    for i=1:length(listing)
        fname = listing(i).name;
        if listing(i).isdir==false
           dot_pos=find(fname=='.');
           if length(dot_pos)>1
               dot_pos=dot_pos(end);
           end
           ext = fname(dot_pos+1:end);
           flag_ext=0;
           if length(ext)==3
               if ext=='csv'
                    flag_ext=1;
               end
           elseif length(ext)==4
               if ext=='xlsx'
                    flag_ext=1;
               end
           end
           if flag_ext==1
               fname_new = strcat(output_folder,fname(1:dot_pos-1),'_res.xlsx')
               opts = detectImportOptions(fname);
%                preview(fname,opts)
               opts.PreserveVariableNames = true;
               DT = readtable(fname,opts);
%                DT.Properties.VariableNames
               DT_new = DT;
               %% remove first 10 lines
               DT_new(1:10,:)=[];
               T = table2array(DT_new(:,1));
               [Tsort, Isort] = sort(T);
               diff_ind=find([1:length(Isort)]'-Isort~=0);
               if ~isempty(diff_ind)
                   warning(sprintf('The order of time instants is disturbed starting from the element #%d',min(diff_ind)+10))
                   DT_new = sortrows(DT_new,'Time [ms]');
                   % find indices of repeating elements
                   idxRepeat = find(diff(Tsort) == 0);
                   DT_new(idxRepeat,:)=[];
               end
               DT_new_L = DT_new(:,[1 5:21]);
               DT_new_R = DT_new(:,[1 22:38]);
%                DT_new_L.Properties.VariableNames
               %% clean blinking
               Pupil_L = table2array(DT_new_L(:,2));
               Pupil_R = table2array(DT_new_R(:,2));
               %
               blink_ext_L=Pupil_L.*[1; Pupil_L(1:end-1)].*[Pupil_L(2:end); 1];
               blink_ext_R=Pupil_R.*[1; Pupil_R(1:end-1)].*[Pupil_R(2:end); 1];
               DT_new_L(blink_ext_L==0,:)=[];
               DT_new_R(blink_ext_R==0,:)=[];
               %% Read and analyse data
               % interpupil distance
               if ~(isempty(DT_new_L)|isempty(DT_new_R))
                   blink_ext = blink_ext_L.*blink_ext_R;
                   DT_new(blink_ext==0,:)=[];
                   InterPupilDist = table2array(DT_new(:,39));
                   [Mean_InterPupilDist,SD_InterPupilDist] = stat(InterPupilDist);
               else
                   Mean_InterPupilDist=0;
                   SD_InterPupilDist=0;
                   if isempty(DT_new_L)
                       DT_new(blink_ext_R==0,:)=[];
                   else 
                       DT_new(blink_ext_L==0,:)=[];
                   end
               end
               % time samples after cleaning
               T_L = table2array(DT_new_L(:,1));
               T_R = table2array(DT_new_R(:,1));
               if ~isempty(T_L)
                   % pupil diameter 
                   PupilDiam_L = table2array(DT_new_L(:,7));
                   %
                   % Average refraction 
                   AR_L = table2array(DT_new_L(:,17));
                   % Gaze: X_L, Y_L
                   X_L = table2array(DT_new_L(:,13));
                   Y_L = table2array(DT_new_L(:,14));
                   %
                   [Mean_L,SD_L,Med_L,Q1_L,Q3_L,IQR_L] = stat(AR_L);
                   [Mean_Pupil_L,SD_Pupil_L] = stat(PupilDiam_L);
                   [Mean_X_L,SD_X_L] = stat(X_L);
                   [Mean_Y_L,SD_Y_L] = stat(Y_L);
               else
                   Mean_L = 0; SD_L = 0; Med_L = 0; Q1_L = 0; Q3_L = 0; 
                   IQR_L = 0; Mean_Pupil_L = 0; SD_Pupil_L = 0;
                   Mean_X_L = 0; SD_X_L = 0; Mean_Y_L = 0; SD_Y_L = 0;
               end
               if ~isempty(T_R)
                   % pupil diameter 
                   PupilDiam_R = table2array(DT_new_R(:,7));
                   %
                   % Average refraction 
                   AR_R = table2array(DT_new_R(:,17));
                   % Gaze: X_R, Y_R
                   X_R = table2array(DT_new_R(:,13));
                   Y_R = table2array(DT_new_R(:,14));
                   %
                   [Mean_R,SD_R,Med_R,Q1_R,Q3_R,IQR_R] = stat(AR_R);
                   [Mean_Pupil_R,SD_Pupil_R] = stat(PupilDiam_R);
                   [Mean_X_R,SD_X_R] = stat(X_R);
                   [Mean_Y_R,SD_Y_R] = stat(Y_R);
               else
                   Mean_R = 0; SD_R = 0; Med_R = 0; Q1_R = 0; Q3_R = 0; 
                   IQR_R = 0; Mean_Pupil_R = 0; SD_Pupil_R = 0;
                   Mean_X_R = 0; SD_X_R = 0; Mean_Y_R = 0; SD_Y_R = 0;
               end
 
               x = cell(1,1); 
               Tbl=table(Mean_L,SD_L,x,Mean_R,SD_R,x,Med_L,Q1_L,Q3_L,IQR_L,x,Med_R,Q1_R,Q3_R,IQR_R,x,...
                   Mean_Pupil_L,SD_Pupil_L,x,Mean_Pupil_R,SD_Pupil_R,x,Mean_InterPupilDist,SD_InterPupilDist,...
                   x,Mean_X_L,SD_X_L,x,Mean_Y_L,SD_Y_L,x,Mean_X_R,SD_X_R,x,Mean_Y_R,SD_Y_R);
               writetable(DT,fname_new,'Sheet','Original data')
               writetable(DT_new,fname_new,'Sheet','Cleaned')
               writetable(Tbl,fname_new,'Sheet','Numerical characteristics')

                %% Fourier analysis
                f_max =2.5;
                if ~isempty(T_L)
                    if iscell(AR_L)
                        AR_L = str2double(AR_L);
                    end
                    AR_Lxx = AR_L-mean(AR_L);
                    Tx = (T_L-T_L(1))*1e-3;
                    %
                    % Left eye
                    [pxx_L,f_L] = plomb(AR_Lxx,Tx,f_max);
                    lf_ind = find(f_L<=0.6);
                    mf_ind = find((f_L>=0.7)&(f_L<=0.9));
                    hf_ind = find((f_L>=1)&(f_L<=2.1));
                    lf_left_mean = mean(pxx_L(lf_ind));
                    lf_left_sd = std(pxx_L(lf_ind));
                    mf_left_mean = mean(pxx_L(mf_ind));
                    mf_left_sd = std(pxx_L(mf_ind));
                    hf_left_mean = mean(pxx_L(hf_ind));
                    hf_left_sd = std(pxx_L(hf_ind));
                else
                    lf_left_mean = 0; lf_left_sd = 0;
                    mf_left_mean = 0; mf_left_sd = 0;
                    hf_left_mean = 0; hf_left_sd = 0;
                end
                if ~isempty(T_R)
                    if iscell(AR_R)
                        AR_R = str2double(AR_R);
                    end
                    AR_Rxx = AR_R-mean(AR_R);
                    Tx = (T_R-T_R(1))*1e-3;
                    % Right eye
                    [pxx_R,f_R] = plomb(AR_Rxx,Tx,f_max);
                    lf_ind = find(f_R<=0.6);
                    mf_ind = find((f_R>=0.7)&(f_R<=0.9));
                    hf_ind = find((f_R>=1)&(f_R<=2.1));
                    lf_right_mean = mean(pxx_R(lf_ind));
                    lf_right_sd = std(pxx_R(lf_ind));
                    mf_right_mean = mean(pxx_R(mf_ind));
                    mf_right_sd = std(pxx_R(mf_ind));
                    hf_right_mean = mean(pxx_R(hf_ind));
                    hf_right_sd = std(pxx_R(hf_ind));
                else
                    lf_right_mean = 0; lf_right_sd = 0;
                    mf_right_mean = 0; mf_right_sd = 0;
                    hf_right_mean = 0; hf_right_sd = 0;
                end
                Tbl_f=table(lf_left_mean, lf_left_sd, mf_left_mean,...
                    mf_left_sd, hf_left_mean, hf_left_sd, x,...
                    lf_right_mean, lf_right_sd, mf_right_mean,...
                    mf_right_sd, hf_right_mean, hf_right_sd);
               writetable(Tbl_f,fname_new,'Sheet','Fourier') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 figure()
%                 plomb(AR_Lx,Tx)
%                 figure()
%                 plot(f,pxx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 % fft
%                 L = 2^(ceil(log2(length(AR_Lxx)))+2); %length of the zero-padded block
%                 AR_L_pad = zeros(1,L);
%                 AR_L_pad(1:length(AR_Lxx)) = AR_Lxx;
%                 res = Fs/L; % resolution
%                 num_el = floor(2.5/res); % number of computed frequencies
%                 Y = fft(AR_L_pad);
%                 P1 = abs(Y(1:num_el)/L);
%                 P1(2:end) = 2*P1(2:end);
%                 ff=[1:length(P1)]*res;
%                 hold on
%                 plot(ff,P1)
           end
        end
    end
    
    %% FUNCTIONS
    function [Mean,SD,Med,Q1,Q3,IQR] = stat(Data)
        if iscell(Data)
            Data=str2double(Data);
        end
        Mean=mean(Data);
        SD=std(Data);
        Q = quantile(Data,[0.25 0.5 0.75]);
        Q1=Q(1);    Med=Q(2);   Q3=Q(3);
        IQR = Q3-Q1;
    end
end