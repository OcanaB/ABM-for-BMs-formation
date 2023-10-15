%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function identify neighbors when a cell is in the border
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function neighbors=borders(grid,i,j,N)
    row0=[i-1 i-1 i-1  i  i  i   i+1 i+1 i+1];
    col0=[j-1  j  j+1 j-1 j j+1  j-1  j  j+1];
    row=row0(and(and(row0~=0,row0~=N+1),and(col0~=0,col0~=N+1)));
    col=col0(and(and(row0~=0,row0~=N+1),and(col0~=0,col0~=N+1)));
    for i=1:length(row)
        neighbors(i)=[grid(row(i),col(i))];
    end
end