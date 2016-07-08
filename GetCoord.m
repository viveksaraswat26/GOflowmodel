function [ a,b,c ] = GetCoord( i,j,totNodes, secNodes )
%Returns an array of neighbour coordinates. First element is left
%node's coordinate, second element is right node's coordinate and third
%element is third node's coordinate
% i,j are the coordinates of node passed, totNodes and secNodes are 

if(ismember(totNodes(i,j),secNodes(i,:)))
    %some code here if its a secondary node
%     if(j>1)
       a = [i,j-1];
%     else
%         a = 0;
%     end
    b = [i,j+1];
    
    for m=1:length(totNodes)
        
        if (totNodes(i-1,m) == totNodes(i,j))
            z=m;
            break;
        end
        
    end
    
    c = [i-1,z];   
    
else
    %For primary node
%     if(j>1)
        a = [i,j-1];
%     else
%         a = 0;
%     end
    b = [i,j+1];
    
    for m=1:length(totNodes)
        
        if (totNodes(i+1,m) == totNodes(i,j))
            z=m;
            break;
        end
        
    end
     
    c = [i+1,z];
    
end



    

end

