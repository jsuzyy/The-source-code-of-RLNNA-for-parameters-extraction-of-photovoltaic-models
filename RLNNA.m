function [BestCost,BestValue,XTarget]=RLNNA(CostFunction,nPop,nVar,VarMin,VarMax,MaxIt,X)


% CostFunction:Objective function which you wish to minimize or maximize
% VarMin:Lower bound of a problem
% VarMax:Upper bound of a problem
% nVar:Number of design variables
% nPop:Population size
% MaxIt:Maximum number of iterations

% BestCost:Convergenc curve
% BestValue: The optimal objective function value
% XTarget:The optimal solution

%% --------------------Initialization----------------------------------------
X_LB=repmat(VarMin,nPop,1);
X_UB=repmat(VarMax,nPop,1);
MaxIt=MaxIt/(2*nPop);
beta=ones(1,nPop); %Intilize modification factor
x_pattern=zeros(nPop,nVar);
cost=zeros(nPop,1);

for i=1:nPop
    x_pattern(i,:)=X(i,:);
    cost(i)=CostFunction(x_pattern(i,:)); %function evaluation
end
[COST,index]=min(cost);
% Creat random initial weights with constraint of Summation each column = 1
ww=ones(1,nPop)*0.5;
w=diag(ww);
for i=1:nPop
    t=rand(1,nPop-1)*0.5;
    t=(t./sum(t))*0.5;
    w(w(:,i)==0,i)=t;
end
%--------------------------------------------------------------------------
XTarget=x_pattern(index,:);   % Best obtained solution
Target=COST;                  % Best obtained objetive function value
wtarget=w(:,index);           % Best obtained weight (weight target)
%% -------------------- Main Loop for  RLNNA -------------------------------
FMIN=zeros(MaxIt,1);
old=x_pattern; %Initilize historical population
BestCost(1)=COST;
for ii=2:MaxIt % Start task
    flag=zeros(1,nPop); 
    if rand<rand % Eq.(32)
        old=x_pattern;
    end
    old=old(randperm(nPop),:); % Eq.(33)
    x_patternxx=x_pattern;
    costxx=cost;
    %------------------ Creating new solutions ----------------------------
    x_new=w*x_pattern;
    x_pattern=x_new+x_pattern;
    %------------------- Updating the weights -----------------------------
    for i=1:nPop
        
        w(:,i)=abs(w(:,i)+((wtarget-w(:,i))*2.*rand(nPop,1))); % Eq.(20)
    end
    
    for i=1:nPop
        w(:,i)=w(:,i)./sum(w(:,i));% Summation of each column = 1
    end
    
    %----------------------- Creat new input solutions by bias operator or transfer operator--------------------
    for i=1:nPop
        
        if rand<beta(i)
            
            %------------- Bias for input solutions -----------------------
            N_Rotate=ceil(beta(i)*nVar);
            
            xx=VarMin+(VarMax-VarMin).*rand(1,nVar);
            rotate_postion=randperm(nVar);rotate_postion=rotate_postion(1:N_Rotate);
            
            for m=1:N_Rotate
                x_pattern(i,rotate_postion(m))=xx(m);%Eq.(22)
            end
            %---------- Bias for weights ----------------------------------
            N_wRotate=ceil(beta(i)*nPop);
            
            w_new=rand(N_wRotate,nPop);
            rotate_position=randperm(nPop);rotate_position=rotate_position(1:N_wRotate);
            
            for j=1:N_wRotate
                w(rotate_position(j),:)=w_new(j,:);%Eq.(23)
            end
            
            for iii=1:nPop
                w(:,iii)=w(:,iii)./sum(w(:,iii));   % Summation of each column = 1
            end
        else
            
            x_pattern(i,:)=x_pattern(i,:)+abs(randn).*(XTarget-old(i,:))+abs(randn).*(XTarget-x_pattern(i,:));%Eq.(31)
            
        end
    end
    x_pattern=max(x_pattern,X_LB);    x_pattern=min(x_pattern,X_UB);     % Check the side constraints
    for i=1:nPop
        cost(i)=CostFunction(x_pattern(i,:));%Function evaluation
        if costxx(i)<cost(i)% Eq.(34)
            cost(i)=costxx(i);
            x_pattern(i,:)=x_patternxx(i,:);
            flag(i)=flag(i)-1;
        end
    end
  %----------------------- Creat new input solutions by feedback operator--------------------  
    for i=1:nPop
        
        a=randperm(nPop,1);
        while a==i
            a=randperm(nPop,1);
        end
        
        
        if cost(a)<cost(i)%Eq.(35)
            X(i,:)=x_pattern(i,:)+abs(randn).*(x_pattern(a,:)-x_pattern(i,:))+abs(randn).*(old(i,:)-x_pattern(i,:));%+gi*rand(1,nVar).*(x_pattern(ind,:)-x_pattern(i,:));
        else
            X(i,:)=x_pattern(i,:)+abs(randn).*(x_pattern(i,:)-x_pattern(a,:))+abs(randn).*(old(i,:)-x_pattern(i,:));%-x_pattern(i,:));
            
        end
        X(i,:) = max(X(i,:), VarMin);
        X(i,:) = min(X(i,:), VarMax);
       tp= CostFunction(X(i,:));
        if tp<cost(i)% Eq.(36)
            cost(i)=tp;
            x_pattern(i,:)=X(i,:);
            flag(i)=flag(i)-1;
        end
        beta(i)=beta(i)*(1+0.01*(flag(i)));% Eq.(27)
        if beta(i)>=1
            beta(i)=1;
        end
    end
    %% ------ Selection ---------------------------------------------------
    
    [FF,Index]=min(cost);
    
    if FF<Target
        Target=FF;
        XTarget=x_pattern(Index,:);
        wtarget=w(:,Index);
    else
        [~,Indexx]=max(cost);
        x_pattern(Indexx,:)=XTarget;
        w(:,Indexx)=wtarget;
    end
    BestCost(ii)=Target;    
end

%% -------------------------------- RLNNA Finishes ----------------------------
BestValue=Target;%Output the optimal solution
end