function [VarMin,VarMax,nVar,X,fun] = PV_Select(func_flag,nPop)
%UNTITLED2 ´Ë´¦ÏÔÊ¾ÓÐ¹Ø´Ëº¯ÊýµÄÕªÒª
%   ´Ë´¦ÏÔÊ¾ÏêÏ¸ËµÃ÷
% From: Wenyin Gong Parameter extraction of solar cell models using repaired adaptive differential evolution

%  PV model with single-diode
if func_flag==1
    VarMin = [ 0.0  0.0      0.0  0.0    1.0];
    VarMax = [ 1.0  1.0e-06      0.5  100.0  2.0];
%     Xmin = [0.760775662  0.323154e-06  0.03637551  53.72563852    1.481225178];
%     Xmax = [0.760775662  0.323154e-06  0.03637551  53.72563852    1.481225178];
    nVar = 5;
    known_optimal   = 0.0;
    for i=1:nPop
        X(i,:)=VarMin+rand(1,nVar).*(VarMax-VarMin);
    end
    fun=@PV_model_single;
end
 
%  PV model with double-diode
if func_flag==2
    VarMin = [ 0.0  0.0      0.0   0.0   1.0  0.0     1.0];
    VarMax = [ 1.0  1.0e-06  0.5  100.0  2.0  1.0e-06  2];
%     Xmin = [0.760781  0.225974e-06  0.03674  55.485441 1.451017 0.749346e-06 2];
%     Xmax = [0.760781  0.225974e-06  0.03674  55.485441 1.451017 0.749346e-06 2];
    nVar = 7;
    for i=1:nPop
        X(i,:)=VarMin+rand(1,nVar).*(VarMax-VarMin);
    end
    fun=@PV_model_double;
end

%  PV model with module-diode
if func_flag==3
    VarMin = [ 0.0  0.0       0.0  0    1.0];
    VarMax = [ 2.0  50.0e-06  2.0  2000.0  50.0];
%     Xmin = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];
%     Xmax = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];    
    nVar = 5;
    for i=1:nPop
        X(i,:)=VarMin+rand(1,nVar).*(VarMax-VarMin);
    end
   fun=@PV_module;
end

if func_flag==4
    VarMin = [ 0.0  0.0       0.0  0    1.0];
    VarMax = [ 2.0  50.0e-06  2.0  2000.0  50.0];
%     Xmin = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];
%     Xmax = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];    
    nVar = 5;
    for i=1:nPop
        X(i,:)=VarMin+rand(1,nVar).*(VarMax-VarMin);
    end
   fun=@PV_module;
end
end

function  result = calculate_objective_single(x,V_L,I_L)

I_ph = x(1);
I_SD = x(2);
R_s	 = x(3);
R_sh = x(4);
n	 = x(5); 
q = 1.60217646e-19;
k = 1.3806503e-23;
T = 273.15 + 33.0;		%  the temperature is set as 33 centi-degree
V_t = k * T / q;
result = I_ph - I_SD * ( exp( (V_L + I_L*R_s) / (V_t*n) ) - 1.0 ) - ( (V_L + I_L*R_s)/R_sh ) - I_L;

end

function obj = PV_model_single(x)

a = load('cell_data.txt');
actual_V_data =  a(:,1);
actual_I_data =  a(:,2);
data_len = length(actual_V_data);
for j=1:data_len
    error_value(j) = calculate_objective_single(x,actual_V_data(j), actual_I_data(j));
end
fitness = sum(error_value.^2);
obj = sqrt(fitness/data_len);

end
% -------------------------  PV model with single-diode
%****************************************************************
% -------------------------  PV module
function  result = calculate_objective_module(x,V_L,I_L)%%%

I_ph = x(1);
I_SD = x(2);
R_s	 = x(3);
R_sh = x(4);
n	 = x(5); 
q = 1.60217646e-19;
k = 1.3806503e-23;
T = 273.15 + 45.0;		%  the temperature is set as 45 centi-degree
V_t = k * T / q;
NS=1;
NP=1;
result = NP*I_ph - NP*I_SD * ( exp( (V_L/NS + (I_L*R_s/NP)) / (V_t*n) ) - 1.0 ) - ( NP*(V_L/NS + (I_L*R_s/NP))/R_sh ) - I_L;

end

function obj = PV_module(x)

a = load('pvmodule_data.txt');
actual_V_data =  a(:,1);
actual_I_data =  a(:,2);
data_len = length(actual_V_data);
for j=1:data_len
    error_value(j) = calculate_objective_module(x,actual_V_data(j), actual_I_data(j));
end
fitness = sum(error_value.^2);
obj = sqrt(fitness/data_len);

end
% -------------------------  PV module
%**********************************************************************
% -----------------------  PV model with double-diode
function  result = calculate_objective_double(x,V_L,I_L)

I_ph	= x(1);
I_SD1	= x(2);
R_s		= x(3);
R_sh	= x(4);
n1		= x(5);
I_SD2	= x(6);
n2		= x(7);
    
q = 1.60217646e-19;
k = 1.3806503e-23;
T = 273.15 + 33.0;		% the temperature is set as 33 centi-degree

result = I_ph - I_SD1 * ( exp( (q*(V_L + I_L*R_s)) / (n1*k*T) ) -1.0 ) - I_SD2 * ( exp( (q*(V_L + I_L*R_s)) / (n2*k*T) ) -1.0 ) - ( (V_L + I_L*R_s)/R_sh ) - I_L;
     
end

function obj = PV_model_double(x)

a = load('cell_data.txt');
actual_V_data =  a(:,1);
actual_I_data =  a(:,2);
data_len = length(actual_V_data);
for j=1:data_len
    error_value(j) = calculate_objective_double(x,actual_V_data(j), actual_I_data(j));
end
fitness = sum(error_value.^2);
obj = sqrt(fitness/data_len);

end
