function angles = genRandomAngle(reduction_factor)

    x = rand(1000,3);
    for i=1:size(x,1)
        x1 = x(i,1);
        x2 = x(i,2);
        x3 = x(i,3);    
        R = rotz(rad2deg(2*pi*x1));
        v = [cos(2*pi*x2)*sqrt(x3); sin(2*pi*x2)*sqrt(x3); sqrt(1-x3)];
        H = eye(3) - 2*(v*v');
        M(:,:,i) = -H*R;
        eulZYX(i,:) = wrapToPi(rotm2eul(M(:,:,i)));    
    end
    
    idx = find(eulZYX(:,2)>(-pi/2/reduction_factor) & eulZYX(:,2)<(pi/2/reduction_factor) ...
             & eulZYX(:,3)>(-pi/2/reduction_factor) & eulZYX(:,3)<(pi/2/reduction_factor  ));
    
    angles = eulZYX(idx(randperm(size(idx,1),1)),:);

end