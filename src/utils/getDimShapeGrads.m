function [Ndim,Bdim] = getDimShapeGrads(N,dN,dim)

switch dim

    case 2
        
        nn2 = 2*length(N);

        Bdim(1,1:2:nn2)  = dN(1,:);
        Bdim(2,2:2:nn2)  = dN(2,:);
        Bdim(3,1:2:nn2)  = dN(2,:);
        Bdim(3,2:2:nn2)  = dN(1,:);

        Ndim(1,1:2:nn2) = N;
        Ndim(2,2:2:nn2) = N;


    case 3
        error('Not yet programmed!')
end

   


end