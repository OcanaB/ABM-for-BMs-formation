%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given a cell (i,j) choosing where to move it or proliferate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,y]=ChooseNeighbor(grid,i,j,N)

%%Exploring neighbors1 of (i,j) cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if or(or(i==1,i==N),or(j==1,j==N))
    row0=[i-1 i-1 i-1  i  i  i   i+1 i+1 i+1];
    col0=[j-1  j  j+1 j-1 j j+1  j-1  j  j+1];
    row=row0(and(and(row0~=0,row0~=N+1),and(col0~=0,col0~=N+1)));
    col=col0(and(and(row0~=0,row0~=N+1),and(col0~=0,col0~=N+1)));
else
    row=[i-1 i-1 i-1  i  i  i   i+1 i+1 i+1];
    col=[j-1  j  j+1 j-1 j j+1  j-1  j  j+1];
end
%ind contains the indices of the neighbors1
ind=sub2ind([N N],row,col);
%neighbors contains the information on what is inside each neighbor positio
if or(or(i==1,i==N),or(j==1,j==N))
    neighbors=borders(grid,i,j,N);
else
    neighbors=[grid(i-1,j-1), grid(i-1,j), grid(i-1,j+1), grid(i,j-1), 1, grid(i,j+1), grid(i+1,j-1), grid(i+1,j), grid(i+1,j+1)];
end
%%next loop allow cells go to the other side of the vessel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:length(neighbors) %For the 8 neighbors of (i,j) (+ itself)
    if neighbors(k)==2 %if it's a vessel
        FreeNeigh=nneigh0(grid,N,ind(k)); %check how many free neighbors this vessel have
        if FreeNeigh==0 %if it has no free neighbors, do nothing
            break
        else %if it has free neighbors
            neighbors(k)=0; %the cell in (i,j) can jump to a free place on the other side
            %choosing where to jump, (a,b) location of the vessel
            [a,b]=ind2sub([N N],ind(k));
            %neighborsAB are neighbors of vessel in (a,b) and indAB 
            neighborsAB=[grid(a-1,b-1), grid(a-1,b), grid(a-1,b+1); grid(a,b-1), 1, grid(a,b+1); grid(a+1,b-1), grid(a+1,b), grid(a+1,b+1)];
            indAB=[ind(k)-N-1, ind(k)-N, ind(k)-N+1, ind(k)-1, ind(k), ind(k)+1, ind(k)+N-1, ind(k)+N, ind(k)+N+1];
            %(a2,b2) contain all positions in neigborsAB which are free
            [a2,b2]=find(neighborsAB==0);
            %choose randomly one of the positions
            r=randi([1 length(a2)],1,1);
            a2=a2(r); b2=b2(r);
            %m is the correspondent indice in indAB
            m=sub2ind([3 3],a2,b2); 
            %Relate indice in andAB with the one in NxN, variable ind
            ind(k)=indAB(m);
        end
    end
end

%%Once no neighbors with vessels are describe, choose where to go
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Nvessel gives how many vessels are around each neighbor of (i,j)
Nvessel=[nneigh2(grid,N,ind(1)) nneigh2(grid,N,ind(2)) nneigh2(grid,N,ind(3)) nneigh2(grid,N,ind(4)) nneigh2(grid,N,ind(5)) nneigh2(grid,N,ind(6)) nneigh2(grid,N,ind(7)) nneigh2(grid,N,ind(8)) nneigh2(grid,N,ind(9)) ];
FreeNeigh=sum(sum(neighbors==0)); %Total number of free neighbors (i,j)
prob=Nvessel.*(neighbors==0); 
FreeNeighWVessel=sum(sum(prob>0)); %# of free neighbors with adjacent vessel
% Probability of going to an space adjacent to the vessel is 80%
prob(prob>0)=0.8/FreeNeighWVessel;
prob(and(neighbors==0,Nvessel==0))=0.2/(FreeNeigh-FreeNeighWVessel);

ind2=datasample(ind,1,'Weights',prob); %choose one of the indices having into account their probabilities
[x,y]=ind2sub([N N],ind2); %gives the position for the new cell/movement
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%nneigh2 give the # of neightbors that are vessels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=nneigh2(grid,N,ind)
    [i,j]=ind2sub([N N],ind);
    if or(or(i==1,i==N),or(j==1,j==N))
        neighbors=borders(grid,i,j,N);
    else
        neighbors=[grid(i-1,j-1), grid(i-1,j), grid(i-1,j+1); grid(i,j-1), 1, grid(i,j+1); grid(i+1,j-1), grid(i+1,j), grid(i+1,j+1)];
    end
    f=sum(sum(neighbors==2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%nneigh0 give the # of free neightbors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=nneigh0(grid,N,ind)
    [i,j]=ind2sub([N N],ind);
    if or(or(i==1,i==N),or(j==1,j==N))
        neighbors=borders(grid,i,j,N);
    else
        neighbors=[grid(i-1,j-1), grid(i-1,j), grid(i-1,j+1); grid(i,j-1), 1, grid(i,j+1); grid(i+1,j-1), grid(i+1,j), grid(i+1,j+1)];
    end
    f=sum(sum(neighbors==0));
end

