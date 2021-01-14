function vec = tosph(poins)
    %vec=zeros(size(poins),size(poins(1,:)));
    vec=zeros(4,3);
    [vec(:,1),vec(:,2),vec(:,3)]=cart2sph(poins(:,1),poins(:,2),poins(:,3));
    vec(:,1)=rad2deg(wrapTo2Pi(vec(:,1)));
    for i=1:size(vec(:,1))
        vec(i,1)=360-vec(i,1);
        if vec(i,1)==360
            vec(i,1)=0;
        end
    end
    vec(:,2)=rad2deg(vec(:,2));
    for i=1:size(vec(:,2))
        if vec(i,2)==90
            vec(i,1)=0;
        end
    end    
end