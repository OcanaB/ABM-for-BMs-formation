%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ABM for the formation of BMs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear

% Initialization
load('vasculature1.mat')
N=80;               % Grid size
grid=vasculature1; 
grid(43,47)=1;grid(42,24)=1;grid(25,36)=1; % Location of initial cells
Ncells=3;           % Number of initial cells
newgrid=grid;
phenot=zeros(N);
newphenot=phenot;

% Parameters
eps=0.035;  % Phenotipic jump
A0=1; A1=5; % Range of Ratio Do-nothig/death
MY=0.05;    % Migration probability

%%
for t=1:50
    if Ncells==0
        disp(['extinct at t=', num2str(t)])
        break
    else
    prolif=0;
    % Identifying where cells are located
    Indexes = find(grid==1);
    Indexes = Indexes(and(and(and(Indexes>2*N,Indexes<N*N-2*N),and(mod(Indexes,N)~=0,mod(Indexes-1,N)~=0)),and(mod(Indexes+1,N)~=0,mod(Indexes-2,N)~=0)));
    % Total number of cells at this time point
    Ncells=length(Indexes);
    % Cell locations in a random order
    Indexes=Indexes(randperm(Ncells));
    % Evaluation of each cell
    for k=1:Ncells
        [i,j] = ind2sub([N N],Indexes(k)); 
        if or(or(i==1,i==N),or(j==1,j==N))
            neighbors=borders(newgrid,i,j,N);
        else
            neighbors=[newgrid(i-1,j-1), newgrid(i-1,j), newgrid(i-1,j+1), newgrid(i,j-1), 1, newgrid(i,j+1), newgrid(i+1,j-1), newgrid(i+1,j), newgrid(i+1,j+1)];
        end
        FreeNeigh=sum(sum(neighbors==0));
        VesselNeigh=sum(sum(neighbors==2));
            % Non-adjacent to vessels death  
            if and(VesselNeigh==0,rand>(phenot(i,j))) 
                newgrid(i,j)=0;
                newphenot(i,j)=0;
            % Free Neighbors    
            elseif FreeNeigh>0      
                random=rand;
                PY=0.35*newphenot(i,j)+0.6;
                A=(A1-A0)*phenot(i,j)+A0;
                betaY=(1-PY-MY)/(A+1);
                % Proliferation
                if random<=PY      
                    [x,y]=ChooseNeighbor(newgrid,i,j,N);
                    newgrid(x,y)=1;
                    newphenot(x,y)=newphenot(i,j)+eps*normrnd(0,1);                    
                    newphenot(i,j)=newphenot(i,j)+eps*normrnd(0,1);
                    prolif=prolif+1;
                % Migration    
                elseif random<=(PY+MY)  
                    [x,y]=ChooseNeighbor(newgrid,i,j,N);
                    newgrid(x,y)=1;
                    newgrid(i,j)=0;
                    newphenot(x,y)=newphenot(i,j)+eps*normrnd(0,1);
                    newphenot(i,j)=0;
                %Death     
                elseif random>(1-betaY)  
                    newgrid(i,j)=0;
                    newphenot(i,j)=0;
                %Do nothing    
                else                    
                    newphenot(i,j)=newphenot(i,j)+eps*normrnd(0,1);
                end
            %No Free Neighbors    
            else
                %Death
                if rand<A/10        
                    newgrid(i,j)=0; 
                    newphenot(i,j)=0;
                %Do nothing    
                else                
                    newphenot(i,j)=newphenot(i,j)+eps*normrnd(0,1);
                end
            end
    end 
    end
    % Variable actualization
    grid=newgrid;
    phenot=newphenot;
    phenot(phenot>1)=1; phenot(phenot<0)=0;
    % Follow-up variables
    FollowGrid(:,:,t)=grid;
    FollowPhenot(:,:,t)=phenot;
    proliferation(t)=prolif;
    Ncells=sum(sum(grid==1));
    % Plot
    imagesc(grid)
    colormap([0 0 0;  48,253,0;  255 0 0]/255);
    axis equal tight off
    set(gca,'YDir','normal')
    drawnow    
    pause(0.01);
end