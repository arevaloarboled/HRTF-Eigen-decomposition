function [ild,itd,inde,ws] = HRinterpol_ws(db,sr,p,s_p)
% db = database with hrtf, need to be in the minimum pahse filter
% sr = samplerate
% p = all positions
% s_p = source_position (AED)
%HRINTERPOL Summary of this function goes here
%   Detailed explanation goes here
load('Delaunay_PKU&IOA_HRTF.mat');
source_pos = round(sphToRect(s_p));
idx = pointLocation(dt,source_pos);
while isnan(idx)
    % dunno why Matlab fails
    tmpP = [source_pos(1)+rand()-0.5,...
            source_pos(2)+rand()-0.5,...
            source_pos(3)+rand()-0.5];
    idx = pointLocation(dt,tmpP);
    %disp('**')
end
c = dt.ConnectivityList(idx,:);
c = dt.Points(c,:);
ws = cartesianToBarycentric(dt,idx,source_pos);
[a,e,d] = cart2sph(c(:,1),c(:,2),c(:,3));
%[a,e,d] = rad_to_deg(a,e,d);
a=wrapTo360(rad2deg(a));
%a=arrayfun(@(x) wrapTo360(x),a);
%a=rad2deg(a);
e=rad2deg(e);
aed = [a,e,d];
inde=[];
for l=1:4
    inde=[inde,find(p(:,3)==round(a(l)) & p(:,2)==round(e(l)) & p(:,1)==round(d(l)))];
end
if sum(ws>1) || isnan(sum(ws))
    ws = ones(1,length(ws))/length(ws);
    %disp('**');
end
h = extractHRIRs(aed);
new_h=[];%zeros(size(h,1),size(h,2),size(h,3));
for i=1:size(h,2)
    for j=1:size(h,3)        
        tmp = resample(h(:,i,j),48000,65536);
        if i==1 && j==1
            new_h=zeros(size(tmp,1),size(h,2),size(h,3));          
        end
        new_h(:,i,j)=tmp;
    end
end
h=new_h;
itds = zeros(4,3);
sr=2^16;
%sr=48000;
for i=1:4
    [itds(i,1),itds(i,2),itds(i,3)] = computeITD(h(:,i,:),sr);
    itds(i,:) = round(itds(i,:)*sr);
end
%mPhase = minPhaseHRIR(h);
%mPhase=[];
ild=[];
for i=1:size(aed,1)
    %ind(p(:,3)==round(aed(i,1)) & p(:,2)==round(aed(i,2)) & p(:,1)==round(aed(i,3)))
    if i==1
        ild=db(:,find(p(:,3)==round(aed(i,1)) & p(:,2)==round(aed(i,2)) & p(:,1)==round(aed(i,3))))*ws(i);
    else
        ild=ild+db(:,find(p(:,3)==round(aed(i,1)) & p(:,2)==round(aed(i,2)) & p(:,1)==round(aed(i,3))))*ws(i);
    end
end
%mPhase = db(find(p(:,3)==a & p(:,2)==e & p(:,1)==d));
%ild = interpolateMinPhase(mPhase,ws);
itd = round(interpolateITDs(itds,ws));

% convolve with sound
%y = [conv(y,h(:,1)), conv(y,h(:,2))];

% apply ITD
%s = [[zeros(t(2),1); y(:,1); zeros(t(3),1)], ...
%     [zeros(t(3),1); y(:,2); zeros(t(2),1)]];

end

