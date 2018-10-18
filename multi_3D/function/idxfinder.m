function [ id1, id2, id3 ] = idxfinder( id, n1, n2, n3 )
        if(id > n1*n2*n3) 
            id = mod(id, n1*n2*n3);
        end
        id1 = mod(id,n1);
        id1(id1 == 0) = n1;
        
        id2 = mod(ceil(id/n1), n2);
        id2(id2 == 0) = n2;
        id3 = ceil(id/(n1 * n2));
end