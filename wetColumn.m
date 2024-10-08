clear variables
%data
IN_DBT=33;  % C
IN_WBT=25; % C
OUT_DBT=[23 24 24.5 26]; % C
OUT_WBT=[18 19 20 20]; % C
Qair=[5 10 15 20] ; % lpm
Qwater= [5 10 20 30]; %lpm
T_in=30; %C
L=1; % m
D=0.025; % m
Area=pi*D^2/4; % m^2
v=Qair.*1.667*10^-5/Area; % m/s
u=1.81*10^-5;% Pa-s
Rho=0.876; % Kg/m^3
Dab=2.49*10^-5; %m^2/s 

%Antoine Parameters
A=8.0713;B=1730.63;C=233.426; 

%RelativeHumidity
IN_RH=RelativeHumidity(IN_DBT,IN_WBT);
OUT_RH=RelativeHumidity(OUT_DBT,OUT_WBT);

%VapourPressure
VP_in=(10.^(A-(B./(C+IN_DBT)))).*133.332; %Pa
VP_out=10.^(A-(B./(C+OUT_DBT))).*133.332; %Pa

%PartialPressure
PartialP_IN=VP_in.*IN_RH; %Pa
PartialP_OUT=VP_out.*OUT_RH; %Pa

%Concentrations
C_IN=PartialP_IN./(8.314.*(IN_DBT+273.15)); %mol/m^3
Wf_in=C_IN*18/1300;
C_OUT=PartialP_OUT./(8.314*(OUT_DBT+273.15)); %mol/m^3
Wf_Out=C_OUT.*18/1300;

%enthapy Calculation

H_in=(1.006*IN_WBT+Wf_in*(1.86*IN_WBT+2501));
H_out=(1.006.*OUT_WBT+Wf_in*(1.86.*OUT_WBT+2501));

% Dimensionless Numbers

Sc=Schmidt_No(u,Rho,Dab);
Re=Reynolds(Rho,v,D,u);
Sh=Sherwood(Re,Sc);

k=Dab.*Sh/D; %m/s

NA=k.*abs(C_OUT-C_IN); %mol/m^2.s
Rate=NA.*pi*D*L; %mol/s
% temperature calculation considering water and air flow rates

Qwater_m3s = Qwater / 60000;
Qair_m3s = Qair / 60000; %
T_out = T_in - ((H_in - H_out) ./ (Qwater_m3s ./ Qair_m3s));

if any(T_out < 0)
    disp('Flow rate too low')
end



% Display results
disp('Output Temperatures: ')
disp(T_out)
disp('Mass Transfer Rates %mol/s: ')
disp(Rate)



% Plotting
figure;

% Plot Qair vs. T_out
subplot(3,2,1);
plot(Qair, T_out, '-o', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'b');
title('Air Flow Rate vs. Output Temperature');
xlabel('Air Flow Rate (lpm)');
ylabel('Output Temperature (°C)');
grid on;

% Plot Qair vs. Rate
subplot(3,2,2);
plot(Qair, Rate, '-s', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'r');
title('Air Flow Rate vs. Mass Transfer Rate');
xlabel('Air Flow Rate (lpm)');
ylabel('Mass Transfer Rate (mol/s)');
grid on;

% Plot Qair vs. Re
subplot(3,2,3);
plot(Qair, Re, '-^', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'g');
title('Air Flow Rate vs. Reynolds Number');
xlabel('Air Flow Rate (lpm)');
ylabel('Reynolds Number (Re)');
grid on;

% Plot Qwater vs.T_out
subplot(3,2,4);
plot(Qwater,T_out ,'-d', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'm');
title(' water Flow Rate vs TemperatureOut.');
xlabel('Water Flow Rate (lpm)');
ylabel('TemperatureOut C');
grid on;

% Plot Qair vs. Sh
subplot(3,2,5);
plot(Qair, Sh, '-p', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'c');
title('Air Flow Rate vs. Sherwood Number');
xlabel('Air Flow Rate (lpm)');
ylabel('Sherwood Number (Sh)');
grid on;

% Plot Qair vs. k
subplot(3,2,6);
plot(Qair, k, '-h', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'y');
title('Air Flow Rate vs. Mass Transfer Coefficient');
xlabel('Air Flow Rate (lpm)');
ylabel('Mass Transfer Coefficient (k) (m/s)');
grid on;