% Matlab version of 'Apoplasm model', written by Cesar Barbedo Rocha, Julie Granger, and Michael Ngari Mathuri

clear all; close all; clc

% Model parameters
epsilon_equil = 15;
% epsilon_equil = 28;
OH = 3e-8;                                                            % OH = 3e-8 at pH 6.5 
% OH = 1e-6;                                                          % OH = 1e-6 at pH 8.0
Mult = 4;
mu_max = 0.8/(24*60*60);		                                      % 1/sec
PC_NH3 = 1.8e-3;				                                      % UNITS: cm/sec
PC_NH4 = 1.8e-6;
alpha_GS = 1.020;
exp_beta = 35.32;
vol_apo = 4/3*3.14*((6.6^3)-(6.5^3))*1e-15;

% integration parameters (in seconds)
dt = 0.0001;                                                          % initial time-step
tmax = 3.75*1640;                                                     % run the simulation until 2.4 days
tsave = 10*60;                                                        % save output every 10 minutes 

% Prepare for saving
filename = 'model_output_2C.txt';
variables_to_save = {'t','Cell_N_total','NH4_cytoplasm', 'C_NH4_apo',...
                     'C_NH3_apo','NH4_out_total','NH3_out','lnNH4_out',...
                     'd15NH4_out','d15NH4apo','UptakeNH414','EffluxNH314',...
                     'EffluxNH414','NH4diffout14','NH3_diffin14','NH4_diffin14',...
                     'd15NH3cyt','d15NH4cyt','CellsL'};
    
% Initial reservoirs
CellN14 = 1.66*(1-.00367);
CellN15 = 1.66*.00367;
NH4_apo14 = 100*(1-.003693)*vol_apo*602409;
NH4_apo15 = 100*.003693*vol_apo*602409;
NH4_in14 = 8.8e-4*(1-.00367);
NH4_in15 = 8.8e-4*.00367;
NH4_out14 = 100*(1-.003693);
NH4_out15 = 100*(.003693);

Reservoirs = [CellN14, CellN15, NH4_apo14, NH4_apo15, NH4_in14 ...
                                     , NH4_in15, NH4_out14, NH4_out15];
                                
% Time-step the model
% initialize iteration variables
i = 0;
t = 0;
flag = [1,1,1,1];

while t<=tmax
    
    i = i+1;
    t = t + dt; 
    
    % Update time-dependent variables
    % Cell_N_total = CellN14+CellN15+NH4_in14+NH4_in15;
    Cell_N_total = Reservoirs(1)+Reservoirs(2);                       % +Reservoirs(5)+Reservoirs(6);
    CellsL = Cell_N_total/1.66e-6;				                      % Reservoir = 1 L
    % NH4cytoplasm = NH4_in14+NH4_in15;
    NH4cytoplasm = Reservoirs(5)+Reservoirs(6);
    
    NH3cytoplasm = NH4cytoplasm*1e-7/1.76e-5;                         % at pH 7
    NH4_cytoplasm = (NH4cytoplasm)/(CellsL*1.15e-12);
    % C_NH4_apo = (NH4_apo14+NH4_apo15)/vol_apo/602409;
    C_NH4_apo = (Reservoirs(3)+Reservoirs(4))/(vol_apo*CellsL);
    NH3_cytoplasm = NH3cytoplasm/(CellsL*1.15e-12);
    C_NH3_apo = C_NH4_apo*OH/1.76e-5;                                 % OH=3e-8 at pH 6.5 (try at pH 8: 1e-6)
    % NH4_out_total = NH4_out15+NH4_out14;
    NH4_out_total = Reservoirs(7)+Reservoirs(8);
    % NH3_out = (NH4_out14+NH4_out15)*1.58e-6/1.76e-5;                % at pH 8.2 (external reservoir = 1 L)
    NH3_out = (Reservoirs(7)+Reservoirs(8))*1.58e-6/1.76e-5;          % at pH 8.2 (external reservoir = 1 L)
    Delta_NH3cyt = (NH3_cytoplasm - C_NH3_apo);                       % moles/cm^3
    Delta_NH4cyt = (NH4_cytoplasm-C_NH4_apo*exp_beta)/(1-exp_beta);
 
    DeltaC_NH3 = NH3_out-C_NH3_apo;
    DeltaC_NH4 = NH4_out_total-C_NH4_apo;
    Q = (1e-5/1000/6.6e-4)*(4*3.14*(6.6e-4)^2)*CellsL;

    MMuptake = C_NH4_apo/(.05+C_NH4_apo);
    MM_GS = NH4_cytoplasm/(NH4_cytoplasm+10);
    ENH3 = PC_NH3*5.3e-6/1000*CellsL;                                 % cm/day/cm^2
    ENH4 = -3.565*PC_NH4*Delta_NH4cyt*5.3e-6/1000;

    % F15NH4out = NH4_out15/(NH4_out15+NH4_out14);
    F15NH4out = Reservoirs(8)/(Reservoirs(8)+Reservoirs(7));
    % F15NH4_apo = NH4_apo15/(NH4_apo14+NH4_apo15);
    F15NH4_apo = Reservoirs(4)/(Reservoirs(3)+Reservoirs(4)); 
    % F15NH4cyt = NH4_in15/(NH4_in15+NH4_in14);
    F15NH4cyt = Reservoirs(6)/(Reservoirs(5)+Reservoirs(6));
    % d15NH4cyt = (((NH4_in15/NH4_in14)/.00367)-1)*1000;
    d15NH4cyt = (((Reservoirs(6)/Reservoirs(5))/.00367)-1)*1000;
    d15NH3cyt = d15NH4cyt-epsilon_equil;                              % epsilon_equil default 15? (try at 28? for apo pH 6.5 & 8.0)
    R15NH3cyt = (d15NH3cyt/1000+1)*.00367;
    F15NH3cyt = R15NH3cyt/(1+R15NH3cyt);
    % d15NH4_out = (((NH4_out15/NH4_out14)/.00367)-1)*1000;
    d15NH4_out = (((Reservoirs(8)/Reservoirs(7))/.00367)-1)*1000;
    d15NH3_out = d15NH4_out-epsilon_equil;
    R15NH3_out = (d15NH3_out/1000+1)*.00367;
    F15NH3_out = R15NH3_out/(1+R15NH3_out);
    % lnNH4_out = -log(NH4_out14+NH4_out15);
    lnNH4_out = -log(Reservoirs(7)+Reservoirs(8));
    % d15NH4apo = (((NH4_apo15/NH4_apo14)/.00367)-1)*1000;
    d15NH4apo = (((Reservoirs(4)/Reservoirs(3))/.00367)-1)*1000;
  
    % Calculate fluxes
    UptakeNH414 = (MMuptake*Mult)*mu_max*Cell_N_total*(1-F15NH4_apo);
    UptakeNH415 = (Mult*MMuptake)*mu_max*F15NH4_apo*Cell_N_total;
    GS14 = mu_max*Cell_N_total*MM_GS*(1-F15NH4cyt);
    GS15 = mu_max*Cell_N_total*MM_GS*F15NH4cyt/alpha_GS;
    EffluxNH314 = Delta_NH3cyt*ENH3*(1-F15NH3cyt);
    EffluxNH315 = Delta_NH3cyt*ENH3*F15NH3cyt;
    EffluxNH414 = ENH4*CellsL*(1-F15NH4cyt);
    % EffluxNH414 = ENH4*CellsL*(1-F15NH4cyt);
    EffluxNH415 = ENH4*CellsL*F15NH4cyt;
    
    if EffluxNH414>0
       EffluxNH414 = EffluxNH414;
       EffluxNH415 = EffluxNH415;
    else
       EffluxNH414 = 0;
       EffluxNH415 = 0;
    end
 
    % DeltaC_NH4 = NH4_out_total-C_NH4_apo;
    if DeltaC_NH4<0
        NH4diffout14 = -DeltaC_NH4*Q*(1-F15NH4_apo);                          
        NH4diffout15 = -DeltaC_NH4*Q*F15NH4_apo;   
    else
        NH4diffout14 = 0;
        NH4diffout15 = 0;
    end
    
    NH3_diffin14 = Q*DeltaC_NH3*(1-F15NH3_out);
    NH3_diffin15 = Q*DeltaC_NH3*F15NH3_out;
  
    if DeltaC_NH4>=0
        NH4_diffin14 = DeltaC_NH4*Q*(1-F15NH4out);
        NH4_diffin15 = F15NH4out*Q*DeltaC_NH4;
    else 
        NH4_diffin14 = 0;
        NH4_diffin15 = 0;
    end
    
    % Assign RHS array to avoid reapeted operations
    RHS1 = GS14;
    RHS2 = GS15;
    RHS3 = EffluxNH314+EffluxNH414+NH3_diffin14+NH4_diffin14... 
                                               -UptakeNH414-NH4diffout14;
    RHS4 = EffluxNH315+EffluxNH415+NH3_diffin15+NH4_diffin15 ...
                                -UptakeNH415-NH4diffout15;                                                                              
    RHS5 = UptakeNH414-GS14-EffluxNH314-EffluxNH414;
    RHS6 = UptakeNH415-GS15-EffluxNH315-EffluxNH415;
    RHS7 = NH4diffout14-NH3_diffin14-NH4_diffin14;
    RHS8 = NH4diffout15-NH3_diffin15-NH4_diffin15;
    
    Fluxes = [RHS1, RHS2, RHS3, RHS4, RHS5, RHS6, RHS7, RHS8];
  
    % First time-step uses forward Euler   
    if i == 1 
        Reservoirs = Reservoirs+Fluxes*dt;
        Fluxes_old = Fluxes;       
    % Second time-step uses second-order Adams-Bashforth
    elseif i == 2
         Reservoirs = Reservoirs+(dt/3)*(3*Fluxes-Fluxes_old); 
         Fluxes_old_old = Fluxes_old;
         Fluxes_old = Fluxes;  
    % Third time-step on use third-order Adams-Bashforth
    else     
         Reservoirs = Reservoirs+(dt/12)*(23*Fluxes-16*Fluxes_old ...
                                                      + 5*Fluxes_old_old);
         Fluxes_old_old = Fluxes_old;
         Fluxes_old = Fluxes;   
    end

    % terminate the loop if the solutions blows up
    if any(isnan(Reservoirs))
        break
    end
    
     % TO DO: implement adaptive timestepping
     if (t/86400 > .75) & flag(1);
          disp('new dt')
          dt = dt/2.
          flag(1) = 0;
      elseif (t/86400 > 1.5) & flag(2);
          disp('new dt')
          dt = dt/2.
          flag(2) = 0;   
      elseif (t/86400 > 2.5) & flag(3);
          disp('new dt')
          dt = dt/2.
          flag(3) = 0;   
      elseif (t/86400 > 3.) & flag(4);
          disp('new dt')
          dt = dt/2.
          flag(4) = 0;   
      end
  
    % Store relevant variables every tsave           
    if i == 1
        
        for j=1:length(variables_to_save)
            eval(['soutput.' variables_to_save{j} '=' variables_to_save{j} ';']);
        end                   

    elseif ~rem(i,fix(tsave/dt))    

        disp(t/86400)

        for j=1:length(variables_to_save)
            eval(['soutput.' variables_to_save{j} '=' '[soutput.' ...
                variables_to_save{j} ',' variables_to_save{j} '];']);
        end
        
    end  
    

end

% Save model output 
fileID = fopen(filename,'w');

% header
header = [string(variables_to_save{1})];
for i=2:length(variables_to_save)
    header = [header, string(variables_to_save{i})];
end
fprintf(fileID,[repmat('%15s ',1,length(header)) '\n'], header);

% data
output = cell2mat(struct2cell(soutput));
fprintf(fileID,[repmat('%15.10f ',1,length(header)) '\n'],output);

fclose(fileID);
                          
