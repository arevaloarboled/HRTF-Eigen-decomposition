function HRTF_Eigen_spatialization(audio,p,output_name,pinna)
% spatialize binaurally the vector audio so its apparent location 
% will be described by the path p and storage the binaural audio in output_name.
%
% SYNOPSIS: HRTF_Eigen_spatialization(audio,p,output_name,ear)
%
% INPUT audio: an monophonic audio vector sampled at 48 kHz
%       p: a vector with the location path in spherical coordinates sepcifing the millisecond in which
%          the sound source will have the described location, first position is taken as a initial position.
%          (distance, elevation, azimuth, seconds)
%       output_name: file name of output to storage the binaural audio
%       file, by default will be "spatialization.wav"
%       pinna: specfify the size of the pinna between Large 'l' or Small 's', by
%       default is taken the Large pinna
%
% SEE ALSO 
%
% AUTHOR    : Camilo Arevalo
% $DATE     : 12-Jun-2020 $
% $Revision : 1.00 $
% DEVELOPED : 9.8.0.1323502 (R2020a)
% FILENAME  : HRTF_Eigen_spatialization
    if nargin==2
        output_name='spatialization.wav';
    end
    if nargin<=3
        pinna='l';
    end
    load('Eigen_HRTF.mat'); %load the database
    w_s=128; %define window size
    buffer=zeros(1024,2); %audio buffer  
    output=[0,0]; %binaural output audio vector
    win=hann(w_s,'periodic').*0.9;
%     [audio,fs]=audioread(audio);
%     if fs~=48000
%         warning('Audio file is not sampled at 49 kHz, performing sampling change...');
%         audio=resample(audio,48000,fs);
%     end
    audio=[zeros(w_s/2,1);audio;zeros(w_s/2,1)]; %padding zeros
    current_p=p(1,1:3);
    [hrtfs,delays]=get_filter(current_p,pinna,Eigen_HRTF);
    p_i=2; %index position on the path p
    for i=1:w_s/2:size(audio,1)
        if i+w_s>size(audio,1) %cut the audio if is outside of the boundarie of the audio array
            break;
        end
        if size(p,1)>p_i
            if i/48000>=p(p_i,4)
                current_p=p(p_i,1:3);
                p_i=p_i+1;
                [hrtfs,delays]=get_filter(current_p,pinna,Eigen_HRTF);
            end
        end
        in_signal=audio(i:i+w_s-1).*win;    
        in_signal=[in_signal;zeros(512-size(in_signal,2),1)];                        
        fin_signal=fft(in_signal,512);
        reconstruction=real(ifft(fin_signal.*hrtfs(:,1)));        
        buffer(delays(1)+1:delays(1)+512,1)=buffer(delays(1)+1:delays(1)+512,1)+reconstruction;
        reconstruction=real(ifft(fin_signal.*hrtfs(:,2)));        
        buffer(delays(2)+1:delays(2)+512,2)=buffer(delays(2)+1:delays(2)+512,2)+reconstruction;
        output=[[output(:,1);buffer(1:w_s/2,1)],[output(:,2);buffer(1:w_s/2,2)]];        
        buffer=[buffer(w_s/2+1:1024,:);zeros(w_s/2,2)];        
    end
    audiowrite(output_name, output(2:size(output,1),:), 48000);
end

function [hrtfs,delays]=get_filter(pos,pinna,db)
    db_idx=find(db.coor(:,1)==pos(1) & db.coor(:,2)==pos(2) & db.coor(:,3)==pos(3)); %check if the measure is in the database
    hrtfs=zeros(512,2);
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
            delays(1)=round(db.l_delays(idl)*48000); %get the delays for the pinna selected for the left ear
            delays(2)=round(db.l_delays(idr)*48000); %get the delays for the pinna selected for the rigth ear
        else
            hrtfs(:,1)=db.r_mean;
            hrtfs(:,2)=db.r_mean;
            v=db.r_V; %get the v for the pinna selected
            idr=db_idx; %index for the right ear
            idl=[pos(1) pos(2) mod(360-pos(3),360)];
            idl=find(db.coor(:,1)==idl(1) & db.coor(:,2)==idl(2) & db.coor(:,3)==idl(3)); %index for the left ear            
            coe_l=db.r_coe{idl}; %get the coefficients for the pinna selected and for left ear
            coe_r=db.r_coe{idr}; %get the coefficients for the pinna selected and for rigth ear
            delays(1)=round(db.r_delays(idl)*48000); %get the delays for the pinna selected for the left ear
            delays(2)=round(db.r_delays(idr)*48000); %get the delays for the pinna selected for the rigth ear
        end               
    else  %interpolation needed
        [d,e,a]=fetch_position(pos); %fetch the position in the mesh of interpolation calculated offline
        if pinna=='l'
            hrtfs(:,1)=db.l_mean;
            hrtfs(:,2)=db.l_mean;
            v=db.l_V; %get the v for the pinna selected
            idl=find(db.inter_coor(:,1)==d & db.inter_coor(:,2)==e & db.inter_coor(:,3)==a); %find the index of the interpolation weigths for the left ear
            idr=[pos(1) pos(2) mod(360-pos(3),360)];
            [d,e,a]=fetch_position(idr); %fetch the position in the mesh of interpolation calculated offline
            idr=find(db.inter_coor(:,1)==d & db.inter_coor(:,2)==e & db.inter_coor(:,3)==a); %find the index of the interpolation weigths for the rigth ear
            ws_l=arrayfun(@(x) db.l_coe{db.idx(idl,x)}.*db.ws(idl,x),[1:4],'un',0); %caluculate the interpolated coefficients to reconstruct hrtf for the left ear
            ws_r=arrayfun(@(x) db.l_coe{db.idx(idr,x)}.*db.ws(idr,x),[1:4],'un',0); %caluculate the interpolated coefficients to reconstruct hrtf for the right ear
            coe_l=zeros(1,max(cellfun(@(x) size(x,2),ws_l)));
            coe_r=zeros(1,max(cellfun(@(x) size(x,2),ws_r)));
            for j=1:4
                coe_l(1:size(ws_l{j},2))=coe_l(1:size(ws_l{j},2))+ws_l{j};
                coe_r(1:size(ws_r{j},2))=coe_r(1:size(ws_r{j},2))+ws_r{j};
            end
            delays(1)=round(sum(db.l_delays(db.idx(idl,:)).*db.ws(idl,:)')*48000); %interpolate the delays for the pinna selected for the left ear
            delays(2)=round(sum(db.l_delays(db.idx(idr,:)).*db.ws(idr,:)')*48000); %interpolate the delays for the pinna selected for the left ear
        else 
            hrtfs(:,1)=db.r_mean;
            hrtfs(:,2)=db.r_mean;
            v=db.r_V; %get the v for the pinna selected
            idr=find(db.inter_coor(:,1)==d & db.inter_coor(:,2)==e & db.inter_coor(:,3)==a); %find the index of the interpolation weigths for the rigth ear
            idl=[pos(1) pos(2) mod(360-pos(3),360)];
            [d,e,a]=fetch_position(idl); %fetch the position in the mesh of interpolation calculated offline
            idl=find(db.inter_coor(:,1)==d & db.inter_coor(:,2)==e & db.inter_coor(:,3)==a); %find the index of the interpolation weigths for the left ear
            ws_l=arrayfun(@(x) db.r_coe{db.idx(idl,x)}.*db.ws(idl,x),[1:4],'un',0); %caluculate the interpolated coefficients to reconstruct hrtf for the left ear
            ws_r=arrayfun(@(x) db.r_coe{db.idx(idr,x)}.*db.ws(idr,x),[1:4],'un',0); %caluculate the interpolated coefficients to reconstruct hrtf for the right ear
            coe_l=zeros(1,max(cellfun(@(x) size(x,2),ws_l)));
            coe_r=zeros(1,max(cellfun(@(x) size(x,2),ws_r)));
            for j=1:4
                coe_l(1:size(ws_l{j},2))=coe_l(1:size(ws_l{j},2))+ws_l{j};
                coe_r(1:size(ws_r{j},2))=coe_r(1:size(ws_r{j},2))+ws_r{j};
            end
            delays(1)=round(sum(db.r_delays(db.idx(idl,:)).*db.ws(idl,:))*48000); %interpolate the delays for the pinna selected for the left ear
            delays(2)=round(sum(db.r_delays(db.idx(idr,:)).*db.ws(idr,:))*48000); %interpolate the delays for the pinna selected for the left ear
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