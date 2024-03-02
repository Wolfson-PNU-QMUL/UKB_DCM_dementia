function [X,Y,Z] = dashed_cylinder(N, R, r1,r2, nsegments)

%nsegments is the number of solid segments
%there is an equal number of invisible segments in between these solid
%segments
%this creates a "dashed cylinder" equivalent of a dashed line
c=0;
nsegments = nsegments*2;
for i = 1:nsegments

    if i ==1
        nr1 = r1;
    else
        nr1 = nr2; %previous end point
    end

    frac = 1./(nsegments-i+1);
    vector = frac*(r2 - nr1);
    nr2 = nr1 + vector;

    %only incorporate every other segment
    if ~mod(i,2) 
        c = c+1;
        [X(:,:,c), Y(:,:,c), Z(:,:,c)] = cylinder2P(N,R,nr1,nr2);
    end
end

end
