function Eigen_HRTF = CreateEigenHRTF(sr,hrtf_size,max_db,total_size)
% binaural spatialization of the audio vector so its apparent location will be 
% described by the path positions and storage the spatialization result in output_name.
%
% SYNOPSIS: HRTF_Eigen_spatialization(audio,position,output_name,ear)
%
% INPUT sr: Sample rate.
%       hrtf_size: Actual size of the HRIR used to create the HRTF
%                  database.
%       max_db: dB treshold for the reconstruction.
%       total_size: Specify the size of the HRIR to be converted to HRTF.
%                   I.e., the HRIR will be padded with zeros until reaching
%                   the total_size.
% OUTPUT Eigen_HRTF: Eigen HRTF decomposition.
%
% SEE ALSO 
%
% AUTHOR    : Camilo Arevalo
% $DATE     : 03-Jan-2021 $
% $Revision : 1.10 $
% DEVELOPED : 9.8.0.1323502 (R2020a)
% FILENAME  : CreateEigenHRTF
    if nargin==3
        off_set=0;
    else
        off_set=total_size-hrtf_size;
    end
    parfor i=1:3
        %just to check the parfor
    end
    %%if you want to perform the diffuse equalization
    disp(['Performing diffuse equalization with ' num2str(sr) ' sampling rate']);
%     curdir=cd;
%     cd 'diffuse equalization/';
%     HRIRS=perform_diffuse_eq(sr);
%     cd(curdir);
%      load('HRIRS48000DiffuseEqualization.mat') %%only for debug
    load('./diffuse equalization/HRIRS.mat');
    HRIRS{1}=resample(HRIRS{1}',sr,65536)';
    HRIRS{2}=resample(HRIRS{2}',sr,65536)';
    Eigen_HRTF={};
    d=[20 30 40 50 75 100 130 160];
    e=[-40:10:90];
    coor=zeros(8*4,3);
    c=1;
    for i=1:8
        for j=1:14
            if j==11
                A=[0:10:359];
            elseif j==12
                A=[0:15:359];
            elseif j==13
                A=[0:30:359];
            elseif j==14
                A=0;
                coor(c,:)=[d(i) e(j) A];
                c=c+1;
                break;
            else
                A=[0:5:359];
            end
            a=size(A);
            for k=1:a(2)
                coor(c,:)=[d(i) e(j) A(k)];
                c=c+1;
            end
        end
    end
    Eigen_HRTF.coor=single(coor);
    for J=1:2
        if J==1
            disp('Eigen decomposition for the left pinna');
        else
            disp('Eigen decomposition for the rigth pinna');
        end        
        hrirs=HRIRS{J};
        %hrirs=HRIRS{2};        
        original_hrirs=hrirs';
        hrirs=minPhaseHRIR(hrirs');        
        hrirs=[hrirs(1:hrtf_size,:);zeros(off_set,6344)];
        total_size=size(hrirs,1);
        fhrirs=fft(hrirs);
        mean_matrix=mean(hrirs,2);
        fmean_matrix=mean(fhrirs,2);
        hrirs_subs=bsxfun(@minus, hrirs, mean_matrix);
        % calculate covariance
        s=cov(hrirs_subs');
        % obtain eigenvalue & eigenvector
        %%is eigen vectors & values!!!
        [V,D]=eig(s);
        eigval=diag(D);
        % sort eigenvalues in descending order
        eigval= eigval(end:-1:1);
        V = fliplr(V);
        fV=fft(V);
        if J==1
            Eigen_HRTF.l_mean=single(fmean_matrix);
            Eigen_HRTF.l_V=single(fV);
        else
            Eigen_HRTF.s_mean=single(fmean_matrix);
            Eigen_HRTF.s_V=single(fV);
        end
        coes={};
        itds={};
%         reverseStr = '';           
        fprintf('Computing coefficeints:\n');
        fprintf(['\n' repmat('•',1,50) '↓ finishing here\n\n']);
        parfor i=1:6344             
            ild=hrirs(:,i);
            itds{i}=computeITD(original_hrirs(:,i),sr);                        
            hm=ild-mean_matrix;                
            coe=arrayfun(@(x) dot(hm,V(:,x)),[1:total_size]);                                                    
            for m=1:total_size
%                 re_const=cell2mat(arrayfun(@(x) coe(x)*V(:,x),[1:m],'un',0));
%                 re_const=mean_matrix+arrayfun(@(x) sum(re_const(x,:)),[1:hrtf_size])';
                re_const=mean_matrix+(coe(1:m)*V(:,1:m)')';
                sd=sptrl_distor(re_const,ild,total_size,sr);
                if sd<max_db                        
                    break;
                end
%                 psercentDone = m;%100 * m / 750;                  
%                 msg = sprintf([['[' repelem('=',round(c/6344*50)),'>',repelem('•',round((1-c/6344)*50)) ']'] 'Coordinate %d %d %d Checking coefficients: %3.1f'],j,k,l,percentDone); %Don't forget this semicolon
%                 fprintf([reverseStr, msg]);   
%                 reverseStr = repmat(sprintf('\b'), 1, length(msg));                
            end
            if mod(i,round(6344/50))==0
                fprintf('\b#\n');
            end            
            coes{i}=single(coe(1:m));
        end
        itds=cell2mat(itds);       
        if J==1
            Eigen_HRTF.l_coe=coes;
            Eigen_HRTF.l_delays=single(itds');
            Eigen_HRTF.l_V=Eigen_HRTF.l_V(:,1:max(arrayfun(@(x) size(Eigen_HRTF.l_coe{x},2),1:6344)));
        else
            Eigen_HRTF.s_coe=coes;
            Eigen_HRTF.s_delays=single(itds');
            Eigen_HRTF.s_V=single(Eigen_HRTF.s_V(:,1:max(arrayfun(@(x) size(Eigen_HRTF.s_coe{x},2),1:6344))));
        end
    end
    disp('Including interpolations....')
    load('DB_ws_10_2_cos2.mat');
    Eigen_HRTF.idx=int32(DB.idx);
    Eigen_HRTF.ws=single(DB.ws);
    Eigen_HRTF.inter_coor=single(DB.coor);
    Eigen_HRTF.sr=sr;
    Eigen_HRTF.max_db=max_db;
    Eigen_HRTF.hrtf_size=hrtf_size;
    disp('Finished!')
    disp([repelem('~',25) 'Statistics' repelem('~',25)]);
    disp('Left pinna # of coefficients:');
    fprintf('Mean: %f\tMax: %d\tMin: %d\tTotal: %d\n',mean(arrayfun(@(x) size(Eigen_HRTF.l_coe{i},2),1:6344)),max(arrayfun(@(x) size(Eigen_HRTF.l_coe{x},2),1:6344)),min(arrayfun(@(x) size(Eigen_HRTF.l_coe{x},2),1:6344)),sum(arrayfun(@(x) size(Eigen_HRTF.l_coe{i},2),1:6344)));
    tmp1=(sum(arrayfun(@(x) size(Eigen_HRTF.l_coe{i},2),1:6344))*4+hrtf_size*4+hrtf_size*4*max(arrayfun(@(x) size(Eigen_HRTF.l_coe{i},2),1:6344))+6344*4)/1024/1024;
    fprintf('Aprox. size for Left pinna: %f MB\n',tmp1);
    disp('Rigth pinna # of coefficients:');
    fprintf('Mean: %f\tMax: %d\tMin: %d\tTotal: %d\n',mean(arrayfun(@(x) size(Eigen_HRTF.s_coe{i},2),1:6344)),max(arrayfun(@(x) size(Eigen_HRTF.s_coe{x},2),1:6344)),min(arrayfun(@(x) size(Eigen_HRTF.s_coe{x},2),1:6344)),sum(arrayfun(@(x) size(Eigen_HRTF.s_coe{i},2),1:6344)));    
    tmp2=(sum(arrayfun(@(x) size(Eigen_HRTF.s_coe{i},2),1:6344))*4+hrtf_size*4+hrtf_size*4*max(arrayfun(@(x) size(Eigen_HRTF.s_coe{i},2),1:6344))+6344*4)/1024/1024;
    fprintf('Aprox. size for Right pinna: %f MB\n',tmp2);    
    tmp3=(2*prod(size(Eigen_HRTF.ws))+4*prod(size(Eigen_HRTF.ws)))/1024/1024;
    fprintf('Aprox. size for interpolation data base: %f MB\n',tmp3);
    disp([repelem(' ',25) repelem('~',10) repelem(' ',25)])
    fprintf('Aprox. size for the eigen data base: %f MB\n',tmp1+tmp2+tmp3);
    disp([repelem('~',25) repelem('~',10) repelem('~',25)]);
%     coor(find([arrayfun(@(x) size(Eigen_HRTF.l_coe{x},2),1:6344)==arrayfun(@(x) size(NEiegen.l_coe{x},2),1:6344)]==0),:)
end

