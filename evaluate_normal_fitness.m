% From: Wenyin Gong Parameter extraction of solar cell models using repaired adaptive differential evolution

function  obj = evaluate_normal_fitness(x,func_flag)

switch  func_flag
    case 1
        % PV model with single-diode
        %  real variables = 5 ; objectives = 1; constraints = 0     
        obj = PV_model_single(x); 
    
    case 2
        %  PV model with double-diode
        %  real variables = 7 ; objectives = 1; constraints = 0   
        obj = PV_model_double(x);  
    case 3
        %  PV module model
        %  real variables = 5 ; objectives = 1; constraints = 0   
        obj = PV_module(x);        
        

 end

end


% -------------------------  PV model with single-diode
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

