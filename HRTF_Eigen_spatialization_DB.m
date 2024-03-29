function HRTF_Eigen_spatialization_DB(audio,positions,output_name,pinna,DB,gain)
% binaural spatialization of the audio vector so its apparent location will be 
% described by the path positions and storage the spatialization result in output_name.
%
% SYNOPSIS: HRTF_Eigen_spatialization(audio,position,output_name,ear)
%
% INPUT audio: A monophonic audio vector sampled at 48 kHz.
%       positions: A vector with the location path in spherical coordinates specifying 
%                  the second in which the sound source will have the described location, 
%                  the first position is taken as an initial position 
%                  (distance, elevation, azimuth, seconds).
%       output_name: file name of output to storage the binaural audio
%                    file, by default will be "spatialization.wav"
%       pinna: specfify the size of the pinna between Large 'l' or Small 's', by
%              default is taken the Large pinna.
%       DB: Specify the database to use.
%       gain: a gain for the audio, by default is 1.
%
% SEE ALSO 
%
% AUTHOR    : Camilo Arevalo
% $DATE     : 12-Jun-2020 $
% $Revision : 1.20 $
% DEVELOPED : 9.8.0.1323502 (R2020a)
% FILENAME  : HRTF_Eigen_spatialization
    if nargin<=2
        output_name='spatialization.wav';
    end
    if nargin<=3
        pinna='l';
    end
    if nargin<=4
        load('Eigen_HRTF.mat'); %load the database
        DB=Eigen_HRTF;
        sr=48000;    
        w_s=128; %define window size
        total_size=512;
    else        
        sr=DB.sr;
        total_size=size(DB.s_mean,1);
        w_s=round(total_size-DB.hrtf_size);
        bf_size=size(DB.s_mean,1)+round(max(DB.s_delays*sr));
    end
    if nargin<=5
        gain=1;
    end
    buffer=zeros(bf_size,2); %audio buffer  
    output=[0,0]; %binaural output audio vector
    win=hann(w_s,'periodic');
%     [audio,fs]=audioread(audio);
%     if fs~=sr
%         warning('Audio file is not sampled at 49 kHz, performing sampling change...');
%         audio=resample(audio,sr,fs);
%     end
    audio=audio.*0.9;
    audio=[zeros(w_s/2,1);audio;zeros(w_s/2,1)]; %padding zeros
    current_p=positions(1,1:3);
    [hrtfs,delays]=get_filter(current_p,pinna,DB,sr,total_size);
    p_i=2; %index position on the path p
    disp('.•.•spatializing.•.•')
    reverseStr = '';
    for i=1:w_s/2:size(audio,1)
        if i+w_s>size(audio,1) %cut the audio if is outside of the boundarie of the audio array
            break;
        end
        if size(positions,1)>p_i
            if i/sr>=positions(p_i,4)
                current_p=positions(p_i,1:3);
                if current_p(2)<-40
                    current_p=-40;
                end
                if round(current_p(3))==360
                    current_p(3)=0;
                end
                if current_p(1)<20
                    current_p(1)=20;
                end
                p_i=p_i+1;   
                %current_p
                [hrtfs,delays]=get_filter(current_p,pinna,DB,sr,total_size);
            end
        end
        in_signal=audio(i:i+w_s-1).*win;    
        in_signal=[in_signal;zeros(total_size-size(in_signal,2),1)];                        
        fin_signal=fft(in_signal,total_size);
        reconstruction=real(ifft(fin_signal.*hrtfs(:,1)));        
        buffer(delays(1)+1:delays(1)+total_size,1)=buffer(delays(1)+1:delays(1)+total_size,1)+reconstruction;
        reconstruction=real(ifft(fin_signal.*hrtfs(:,2)));        
        buffer(delays(2)+1:delays(2)+total_size,2)=buffer(delays(2)+1:delays(2)+total_size,2)+reconstruction;
        output=[[output(:,1);buffer(1:w_s/2,1)],[output(:,2);buffer(1:w_s/2,2)]];        
        buffer=[buffer(w_s/2+1:bf_size,:);zeros(w_s/2,2)];
        
        msg = sprintf(['[' repelem('=',round(i/size(audio,1)*50)),'>',repelem('•',round((1-i/size(audio,1))*50)) ']']); %Don't forget this semicolon
        fprintf([reverseStr, msg]);   
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    disp('finished!')
    audiowrite(output_name, output(2:size(output,1),:)*gain, sr);
end

function [hrtfs,delays]=get_filter(pos,pinna,db,sr,total_size)
%     sr=48000;
    db_idx=find(db.coor(:,1)==pos(1) & db.coor(:,2)==pos(2) & db.coor(:,3)==pos(3)); %check if the measure is in the database
    hrtfs=zeros(total_size,2);
    delays=zeros(1,2);
    if size(db_idx,1)==1 %no interpolation needed
        if pinna=='l'
            hrtfs(:,1)=db.l_mean;
            hrtfs(:,2)=db.l_mean;
            v=db.l_V; %get the v for the pinna selected
            idl=db_idx; %index for the left ear
            idr=[pos(1) pos(2) mod(360-pos(3),360)];
            idr=find(db.coor(:,1)==idr(1) & db.coor(:,2)==idr(2) & db.coor(:,3)==idr(3)); %index for the right ear                        
            coe_l=db.l_coe{idl}; %get the coefficients for the pinna selected and for left ear
            coe_r=db.l_coe{idr}; %get the coefficients for the pinna selected and for rigth ear
            delays(1)=round(db.l_delays(idl)*sr); %get the delays for the pinna selected for the left ear
            delays(2)=round(db.l_delays(idr)*sr); %get the delays for the pinna selected for the rigth ear
        else
            hrtfs(:,1)=db.r_mean;
            hrtfs(:,2)=db.r_mean;
            v=db.r_V; %get the v for the pinna selected
            idr=db_idx; %index for the right ear
            idl=[pos(1) pos(2) mod(360-pos(3),360)];
            idl=find(db.coor(:,1)==idl(1) & db.coor(:,2)==idl(2) & db.coor(:,3)==idl(3)); %index for the left ear            
            coe_l=db.r_coe{idl}; %get the coefficients for the pinna selected and for left ear
            coe_r=db.r_coe{idr}; %get the coefficients for the pinna selected and for rigth ear
            delays(1)=round(db.r_delays(idl)*sr); %get the delays for the pinna selected for the left ear
            delays(2)=round(db.r_delays(idr)*sr); %get the delays for the pinna selected for the rigth ear
        end               
    else  %interpolation needed
        [d,e,a]=fetch_position(pos); %fetch the position in the mesh of interpolation calculated offline
        if pinna=='l'
            hrtfs(:,1)=db.l_mean;
            hrtfs(:,2)=db.l_mean;
            v=db.l_V; %get the v for the pinna selected
            idl=find(db.inter_coor(:,1)==d & db.inter_coor(:,2)==e & round(db.inter_coor(:,3))==round(a)); %find the index of the interpolation weigths for the left ear
            idr=[pos(1) pos(2) mod(360-pos(3),360)];
            [d,e,a]=fetch_position(idr); %fetch the position in the mesh of interpolation calculated offline
            idr=find(db.inter_coor(:,1)==d & db.inter_coor(:,2)==e & round(db.inter_coor(:,3))==round(a)); %find the index of the interpolation weigths for the rigth ear
            ws_l=arrayfun(@(x) db.l_coe{db.idx(idl,x)}.*db.ws(idl,x),[1:4],'un',0); %caluculate the interpolated coefficients to reconstruct hrtf for the left ear
            ws_r=arrayfun(@(x) db.l_coe{db.idx(idr,x)}.*db.ws(idr,x),[1:4],'un',0); %caluculate the interpolated coefficients to reconstruct hrtf for the right ear
            coe_l=zeros(1,max(cellfun(@(x) size(x,2),ws_l)));
            coe_r=zeros(1,max(cellfun(@(x) size(x,2),ws_r)));
            for j=1:4
                coe_l(1:size(ws_l{j},2))=coe_l(1:size(ws_l{j},2))+ws_l{j};
                coe_r(1:size(ws_r{j},2))=coe_r(1:size(ws_r{j},2))+ws_r{j};
            end
            delays(1)=round(sum(db.l_delays(db.idx(idl,:)).*db.ws(idl,:)')*sr); %interpolate the delays for the pinna selected for the left ear
            delays(2)=round(sum(db.l_delays(db.idx(idr,:)).*db.ws(idr,:)')*sr); %interpolate the delays for the pinna selected for the left ear
        else 
            hrtfs(:,1)=db.r_mean;
            hrtfs(:,2)=db.r_mean;
            v=db.r_V; %get the v for the pinna selected
            idr=find(db.inter_coor(:,1)==d & db.inter_coor(:,2)==e & round(db.inter_coor(:,3))==round(a)); %find the index of the interpolation weigths for the rigth ear
            idl=[pos(1) pos(2) mod(360-pos(3),360)];
            [d,e,a]=fetch_position(idl); %fetch the position in the mesh of interpolation calculated offline
            idl=find(db.inter_coor(:,1)==d & db.inter_coor(:,2)==e & round(db.inter_coor(:,3))==round(a)); %find the index of the interpolation weigths for the left ear
            ws_l=arrayfun(@(x) db.r_coe{db.idx(idl,x)}.*db.ws(idl,x),[1:4],'un',0); %caluculate the interpolated coefficients to reconstruct hrtf for the left ear
            ws_r=arrayfun(@(x) db.r_coe{db.idx(idr,x)}.*db.ws(idr,x),[1:4],'un',0); %caluculate the interpolated coefficients to reconstruct hrtf for the right ear
            coe_l=zeros(1,max(cellfun(@(x) size(x,2),ws_l)));
            coe_r=zeros(1,max(cellfun(@(x) size(x,2),ws_r)));
            for j=1:4
                coe_l(1:size(ws_l{j},2))=coe_l(1:size(ws_l{j},2))+ws_l{j};
                coe_r(1:size(ws_r{j},2))=coe_r(1:size(ws_r{j},2))+ws_r{j};
            end
            delays(1)=round(sum(db.r_delays(db.idx(idl,:)).*db.ws(idl,:)')*sr); %interpolate the delays for the pinna selected for the left ear
            delays(2)=round(sum(db.r_delays(db.idx(idr,:)).*db.ws(idr,:)')*sr); %interpolate the delays for the pinna selected for the left ear
        end        
    end
    %reconstruction    
    hrtfs(:,1)=hrtfs(:,1)+(coe_l*v(:,1:size(coe_l,2))')';
    hrtfs(:,2)=hrtfs(:,2)+(coe_r*v(:,1:size(coe_r,2))')';    
end

function [d,e,a]=fetch_position(pos)
    e=[-40:2:90];
    e=e(round(pos(2)/2)+21);
    d=[20 30 40 50 75 100 130 160];
    for k=1:7
        if pos(1)>=d(k) && pos(1)<=d(k+1) 
            pos(1)=pos(1)-d(k);
            pos(1)=round(pos(1)/((d(k+1)-d(k))/10));
            d=d(k):(d(k+1)-d(k))/10:d(k+1);
            d=d(pos(1)+1);
            break;
        end
    end
    a=360/round(cos(deg2rad(e))*360);
    pos(3)=round(pos(3)/a);    
    a=0:a:359;
    a=a(pos(3)+1);
end