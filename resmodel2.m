% Design a flow model through graphene oxide. Created by: Vivek Saraswat, Materials
% Science and Engineering, University of Wisconsin Madison, June 2016

clear all
close all

h = waitbar(0,'Starting code execution...');

%% create randomly sized sheets from a normal distribution of mean mu and std dev sigma
mu=1; %in micrometers
sigma=mu/5; %in micrometers
nsheets = 1000; %number of sheets per layer
nlayers = 10; %number of layers
sheetLen = normrnd(mu,sigma,nlayers,nsheets); %Create a notmally distributed matrix representing the sheet lengths
GConst = 10; % conductivity per unit lenght of the sheets
toplot = 0;


%Assign coordinates of these sheets. End of a sheet is its coordinate and
%plot it

%% Mark Primary Nodes & identify last nodes
for i=1:nlayers
    
    y=nlayers-i+1;
    
    for j=1:nsheets
        if (j==1)
            primSheetPos(i,j)=sheetLen(i,j);
            if (toplot)
                plot([0 primSheetPos(i,j)],[y y],'bl');
                hold on
            end
            
        else
            primSheetPos(i,j)=primSheetPos(i,j-1)+sheetLen(i,j);
            x1=primSheetPos(i,j-1);
            x2=primSheetPos(i,j);
            if (toplot)
                plot ([x1 x2],[y y],'bl');
            end
                                    
        end
        
        if(toplot)
            plot(primSheetPos(i,j),y,'b+');
        end
               
    end
    
    
    lastNode(i)=primSheetPos(i,nsheets);
    
end


%% Mark Secondary Nodes

for i=2:nlayers
    secSheetPos (i,:) = primSheetPos (i-1,:);
    
    if (toplot)
        plot (secSheetPos(i,:),nlayers-i+1,'r--o');
    end
end

totSheetPos = [primSheetPos secSheetPos];

for i=1:nlayers
    totSheetPos(i,:) = sort(totSheetPos(i,:));
end

%% Matrix Modification
% Modify the totSheetPos matrix to account for the first row
for i=1:nlayers-1
    newtotSheetPos(i,:)=totSheetPos(i+1,:);    
end

% Modify secondary node matrix to remove zeros from layer 1
newSecPos = secSheetPos((2:end),:);

%% Conductivity matrix

n=length(newtotSheetPos);

G = zeros(nlayers-1,n);

for i = 1:nlayers-1
    for j=1:n
        if(j==1)
            G(i,j)=0;%1e-5;%setting conductivity of first sheet as ~ 0 because no current flows out
        else
            var = newtotSheetPos(i,j)- newtotSheetPos(i,j-1);
            G(i,j) = GConst*(1/var)*(newtotSheetPos(i,j) <= lastNode(i+1));%1e-5; % Conductivity of all nodes > last primary node is set to zero
        end
    end
end

%% Writing kirchoff current law at (length(newtotSheetPos)-1)*(nlayers-1) nodes

waitbar(0,h,'Starting to assign the linear equations...')

inpVolt=5; %Choose the input voltage
Gint = GConst/0.001; %Interlayer Conductivity/um / thickness
m=0;
l=1;

coeff = zeros(nlayers-1,n);

for i=1:nlayers-1
    
    for j=1:n
        
        m=m+1;
        
        if(i==1) %First Layer
            
            if(j==1) % First node
                
                if(ismember(newtotSheetPos(i,j),newSecPos(i,:))) %If first node is a secondary node
                    coeff(m,(i-1)*n+j)=1;
                    coeff(m,(i-1)*n+j+1)=-G(i,j+1)/(G(i,j+1)+Gint+G(i,j));
                    b(m)=(Gint+G(i,j))*inpVolt/(G(i,j+1)+Gint+G(i,j)); %Gleft and Gtop are constants...
                
                else
                    coeff(m,(i-1)*n+j)=1;
                    coeff(m,(i-1)*n+j+1)=-G(i,j+1)/(G(i,j+1)+G(i,j));
                    b(m)=G(i,j)*inpVolt/(G(i,j+1)+G(i,j));
                    
                end
                
                
            elseif (j<n)
                
                if(ismember(newtotSheetPos(i,j),newSecPos(i,:)))
                
                    coeff(m,(i-1)*n+j)=1;
                    coeff(m,(i-1)*n+j+1)=-G(i,j+1)/(G(i,j+1)+Gint+G(i,j-1));
                    coeff(m,(i-1)*n+j-1)=-G(i,j-1)/(G(i,j+1)+Gint+G(i,j-1));
                    b(m)=Gint*inpVolt/(Gint+G(i,j+1)+G(i,j-1));
                    
                else
                    [p q r]= GetCoord(i,j,newtotSheetPos,newSecPos);
                    coeff(m,(i-1)*n+j)=1;
                    coeff(m,(i-1)*n+j+1)=-G(i,j+1)/(G(i,j+1)+Gint+G(i,j-1));
                    coeff(m,(i-1)*n+j-1)=-G(i,j-1)/(G(i,j+1)+Gint+G(i,j-1));
                    coeff(m,(r(1)-1)*n+r(2))=-Gint/(G(i,j+1)+Gint+G(i,j-1));
                    b(m)=0;
                    
                end
                
                
            else
                
                coeff(m,(i-1)*n+j)=1;
                b(m)=inpVolt;%/100;
                
            end
        
            
        elseif (i~=nlayers-1)
            
            if(j==1)
                
                [p q r]= GetCoord(i,j,newtotSheetPos,newSecPos);
                coeff(m,(i-1)*n+j)=1;
                coeff(m,(i-1)*n+j+1)=-G(i,j+1)/(G(i,j+1)+Gint+G(i,j));
                coeff(m,(r(1)-1)*n+r(2))=-Gint/(G(i,j+1)+Gint+G(i,j));
                b(m)=G(i,j)*inpVolt/(G(i,j+1)+Gint+G(i,j));
                
            elseif (j<n)
                
                [p q r]= GetCoord(i,j,newtotSheetPos,newSecPos);
                coeff(m,(i-1)*n+j)=1;
                coeff(m,(i-1)*n+j+1)=-G(i,j+1)/(G(i,j+1)+Gint+G(i,j));
                %new addition
                coeff(m,(i-1)*n+j-1)=-G(i,j)/(G(i,j+1)+Gint+G(i,j));
                coeff(m,(r(1)-1)*n+r(2))=-Gint/(G(i,j+1)+Gint+G(i,j));
                b(m)=0;
                
            else
                
                coeff(m,(i-1)*n+j)=1;
                b(m)=inpVolt/100;
                
            end  
            
        elseif (i==nlayers-1) %last layer
            
            if(ismember(newtotSheetPos(i,j),newSecPos(i,:)))
                
                if(j~=1 && j~=n) %* if the node is a secondary node and is not a starting node
                
                    [p q r]= GetCoord(i,j,newtotSheetPos,newSecPos);
                    coeff(m,(i-1)*n+j)=1;
                    coeff(m,(r(1)-1)*n+r(2))=-Gint/(G(i,j+1)+Gint+G(i,j));
                    b(m)=0;
                    
                elseif (j==1)
                    
                    [p q r]= GetCoord(i,j,newtotSheetPos,newSecPos);
                    coeff(m,(i-1)*n+j)=1;
                    coeff(m,(r(1)-1)*n+r(2))=-Gint/(G(i,j+1)+Gint+G(i,j));
                    b(m)=G(i,j)*inpVolt/(G(i,j+1)+Gint+G(i,j));;
                    
                else
                    
                    coeff(m,(i-1)*n+j)=1;
                    b(m)=inpVolt/100;
                    
                    
                end
                
                
            else
                
                
                coeff(m,(i-1)*n+j)=1;
                b(m)=0;
                      
              
                
            end
         
        end
      
    end
    
    waitbar(i/(nlayers-1),h,sprintf('%.2f%% of layers done...',i*100/(nlayers-1)))
end

B=b.';
waitbar(0,h,'now solving linear eqns')
X=linsolve(coeff,B);
Xmat=reshape(X,[ n ,nlayers-1]);
Xmat=Xmat.';


sum=0;

for i=1:n
    if(ismember(newtotSheetPos(1,i),newSecPos(1,:)))
    
        sum=sum + Gint*(inpVolt-X(i));
        
    end
    
end

waitbar(1,h,'Success!')

Req=inpVolt/sum;

disp(Req)

close(h)
% disp(Xmat)

% Req(z)=inpVolt/sum;
% 
% waitbar(z/zmax,h)



