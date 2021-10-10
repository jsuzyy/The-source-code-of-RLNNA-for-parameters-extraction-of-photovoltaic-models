clc
clear all

Number=3;% the number of models
LoopIter=30;% the number of runs
nPop=40;%population size
for index=1:Number % select model: index=1, single diode model; index=2, double diode mdoel; index=3, PV module model
    if index==1
        disp(['Perform NNA and RLNNA for extracting the parameters of single diode model']);
    elseif index==2
        disp(['Perform NNA and RLNNA for extracting the parameters of double diode model']);
    else
        disp(['Perform NNA and RLNNA for extracting the parameters of PV module model']);
    end
    for i=1:LoopIter % start task
        [VarMin,VarMax,nVar,X,fun] = PV_Select(index,nPop);%intilize parameters
        if index==1
            MaxIt=20000;%the maximum number of function evaluations
            [BestCost1_1(i,:),BestValue1_1(i),Best1_1(i,:)]=NNA(fun,nPop,nVar,VarMin,VarMax,MaxIt,X);
            disp(['Iteration ' num2str(i) ': the optimal solution of NNA   = ' num2str(BestValue1_1(i),15)]);
            [BestCost2_1(i,:),BestValue2_1(i),Best2_1(i,:)]=RLNNA(fun,nPop,nVar,VarMin,VarMax,MaxIt,X);
            disp(['Iteration ' num2str(i) ': the optimal solution of RLNNA = ' num2str(BestValue2_1(i),15)]);
        elseif index==2
            MaxIt=50000;%the maximum number of function evaluations
            [BestCost1_2(i,:),BestValue1_2(i),Best1_2(i,:)]=NNA(fun,nPop,nVar,VarMin,VarMax,MaxIt,X);
            disp(['Iteration ' num2str(i) ': the optimal solution of NNA   = ' num2str(BestValue1_1(i),15)]);
            [BestCost2_2(i,:),BestValue2_2(i),Best2_2(i,:)]=RLNNA(fun,nPop,nVar,VarMin,VarMax,MaxIt,X);
            disp(['Iteration ' num2str(i) ': the optimal solution of RLNNA = ' num2str(BestValue2_1(i),15)]);
        else
            MaxIt=22000;%the maximum number of function evaluations
            [BestCost1_3(i,:),BestValue1_3(i),Best1_3(i,:)]=NNA(fun,nPop,nVar,VarMin,VarMax,MaxIt,X);
            disp(['Iteration ' num2str(i) ': the optimal solution of NNA   = ' num2str(BestValue1_1(i),15)]);
            [BestCost2_3(i,:),BestValue2_3(i),Best2_3(i,:)]=RLNNA(fun,nPop,nVar,VarMin,VarMax,MaxIt,X);
            disp(['Iteration ' num2str(i) ': the optimal solution of RLNNA = ' num2str(BestValue2_1(i),15)]);
        end
    end
end