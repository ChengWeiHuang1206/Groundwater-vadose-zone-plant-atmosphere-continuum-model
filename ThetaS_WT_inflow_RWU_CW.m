clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Date: 9/02/2021
% Author: Cheng-Wei Huang
% Goal: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Setting grids for time and space domain and soil-root-plant properties
%---- Computational Grid (Time = second, Length = m)
% Duration:
day=15;
duration=day*24*60*60; %[s]
dt=0.02; %[s]
Tf=duration/dt;  % # of time steps

N_30=30*60/dt; % # of data for 30 mins
N_10=10*60/dt; % # of data for 10 mins
N=Tf/N_30;  % # of 30mins for whole duration
% Tf=186492;  % 186495 turning point

% Tf=360000;
%--------------- Soil------------------
% soil-root properties
soil=11;
%---- Hydraulic Parameters from Clapp and Hornberger (1978)
[b,thetas, ks, psis]=soil_hydraulic_values_generator(soil); %ks in m/s, psis in m

% soil grid
As=1;  %Ground area [m^2]
ms=1600;  % Number of grid
zmins=0;
zmaxs=1.2; %Depth of soil domain for "calculation" [m]
dzs=(zmaxs-zmins)/ms; % Depth of grid [m]
zs=[1:ms]*dzs-dzs/2; % grids: 1 3 5 7 9 ........
% Computational Grid
m=ms;
dz=dzs;
z=zs;
% outer water level
zout=zs(1000);  % [m]
% Depth of soil domain
Ls=5; %[m]
%--------------- Roots ------------------
Lr=1.2; % Rooting depth [m]
Nr=length(find(zs<=Lr)); % Number of vertical grids within root zone
% Root length density: root length per unit soil volume [m m^-3]  {1 X Nr}
lambda=-(0.976.^(zs(1:Nr)*100))*log(0.976)*4000*201.56455685; % Power-law
lambda=[flip(lambda) 0*zs(find(zs>Lr))]; % Reverse power-law

% Root membrane permeability [s^-1]
Kr=10^(-9);
% Effective root radius [m]
rr=3*10^(-3);
% Root surface density: root surface per unit soil volume [m^2 m^-3]  {1 X Nr}
B=2*pi*rr*lambda;
% Length scale characterizing the mean radial distance for the movement of
% water molecules from the bulk soil to the root surface within the rhizosphere
l=0.53./(pi*lambda).^(1/2); %[m]
%--------------- Leaf ------------------
LAI=1.5; %[m^2 m^-2]
Al=LAI*As; % leaf area [m^2]
d=1*10^(-4); % characteristic leaf dimension [m]

% parameters for MWUE(La) (Manzoni et al.,2011)
gcut=0.04;  % Nocturnal residuel conductance
if gcut>0
    La_type=2;
else
    La_type=1;
end
%%  Seting environmental factors
load('Env.mat') % Env=[X_time; Ta; RH; PPFD; U];
Ta=Env(2,:);
RH=Env(3,:);
PPFD=Env(4,:);
U=Env(5,:);
[estar,VPD,ea,rhov,Molar_Conc,VPD_molConc]=vapor_pressure_deficit (Ta,RH);
Ta=repmat(Ta,1,day);
ea=repmat(ea,1,day);
PPFD=repmat(PPFD,1,day);
U=repmat(U,1,day);
VPD=repmat(VPD,1,day);

Pa=101.3; % [kPa]
ca=400; % [ppm]
%%  Setting Initial condition and boundary condition
% Set initial condition for soil water profile in the unsaturated zone
% Analytical solution for the soil water profile of unsaturated zone
% (1) Total head at WT: Matric potential is not zero because Clapp and Hornberger
%                   is used.
z_inner=zs(1000); % Initial groundwater level [m]

z_un=z(max(find(z<=z_inner)));
p_WT=-max(z_un)+psis.*(thetas/thetas).^(-b);
% (2) Total head across the unsaturated zone:
% p_US=-z+psis.*(thetai/thetas).^(-b);
% (1)=(2)
% => p_WT=-z+psis.*(thetai/thetas).^(-b)
% => psis.*(thetai/thetas).^(-b)= p_WT+z
% => (thetai/thetas).^(-b)= (p_WT+z)/psis
% => (thetai/thetas)= ((p_WT+z)/psis).^(-1/b)
thetai =thetas*((p_WT+z(find(z<=z_inner)))/psis).^(-1/b);
thetai(length(thetai))=thetas;

% figure
% subplot(1,2,1)
% plot(thetai,-z(find(z<=z_inner)))
% xlabel('\theta_s (m^{3} m^{-3})')
% ylabel('depth (m)')
% subplot(1,2,2)
% p_ana=-z(find(z<=z_inner))+psis.*(thetai/thetas).^(-b);
% plot(p_ana,-z(find(z<=z_inner)))
% xlabel('head (m)')
% ylabel('depth (m)')

theta_pre=thetai; clear thetai
thetai=ones(1,ms)*thetas;
thetai(1:length(theta_pre))=theta_pre;
thetai_sum=sum(thetai)*dzs;  % Initial total water volumn [m3]
% thetai=smooth(thetai,10)';
% thetass=max(thetai);
% figure
% subplot(1,2,1)
% hold on
% plot(thetai,-z,'k-')
% subplot(1,2,2)
% hold on
% p1=-z+psis.*(thetai/thetas).^(-b); %[m]
% p1(find(thetai==thetas))=-z(min(find(thetai==thetas)))+psis.*(thetas/thetas).^(-b);
% plot(p1,-z,'k-')

%---- Specify upper boundary for soil domain
qs_up=0; %[m s^-1]

%----Set leaf water potential for the past 24 hrs
load('psisL_30avg_oneday.mat')


%%  Space for saving data
TIME=zeros(N,1); % Time [hr]
z_GT=zeros(N,1); % Groundwater table [m]
QQ=zeros(N,1); % Recharge rate: upward flux [m s-1]
psi_L=zeros(N,1); % Leaf water potential [Pa]
psi_B=zeros(N,1); % Water potential at the stem base [Pa]
EXFL=zeros(N,1); % Exfiltration rate

psisL_avg_mm=zeros(N,1); % Averaged leaf water potential for the past 24 hrs [MPa]
La_mm=zeros(N,1); % Marginal water use efficiency
fe_mm=zeros(N,1);
An_mm=zeros(N,1);
gs_mm=zeros(N,1);

theta_mm=zeros(N,ms);
Tr_layer_mm=zeros(N,ms);
Tr_Total_mm=zeros(N,1);
% QQ=zeros(Tf,1);
% z_GT=zeros(Tf,1);
% psi_L=zeros(Tf,1);
% psi_B=zeros(Tf,1);
%%
V1=0;V2=0;

tic
Q_total1=0;
Q_total=0;
F_e_total=0;

for kk=1:N    
    
    kk
    %----- Leaf gas exchange
    % Step 1: Compute leaf-level gas exchange based on stomatal optimization theory using known marginal water
    % use efficiency from previous water status
    psisL_avg_mm(kk)=mean(psisL_30avg_oneday); %[MPa]
    [La_mm(kk)] = MWUE_Manzoni_2011(psisL_avg_mm(kk),ca,Pa,La_type);
    
    if PPFD(kk)==0
        fe_mm(kk)=gcut*VPD(kk)/Pa;
        An_mm(kk)=0;
        gs_mm(kk)=0;
    else
        [cs,gs_mm(kk),gtc,gtv,gb_T,An_mm(kk),fe_mm(kk),H,ci,ei_v,T_s ] = optimized_two_layer_numerical( ca,PPFD(kk),ea(kk),U(kk),d,Ta(kk),Pa, La_mm(kk),gcut);
    end
    
    Exfl=0;
    for i=1:N_30
        
        % Step 2: Determine recharge rate from /discharge rate to external
        % water body and vertical water potentail gradient based on known 
        % soil water status (unsaturated zone) and known water table and depth of external water body 
        %------------------------------------------------------------------
        % Groundwater recharge rate
        N_saturate=min(find(thetai==thetas));
        % Depth of water table [m]
        zw=z(N_saturate);
        % Recharge to the soil [m s-1]
        Q_re=ks*(-zout+zw)./(Ls-zw);
        % Normalized recharge [dimensionless]
        C=Q_re/ks;
        %------------------------------------------------------------------
        % Total water potential in unsaturated zone
        p=-z+psis.*(thetai/thetas).^(-b); %[m]
        % Total water potential in saturated zone     
        p(find(thetai==thetas))=-z(find(thetai==thetas))+(1+C).*z(find(thetai==thetas))-zout-C*Ls+psis; %[m]            
        %------------------------------------------------------------------
        % Soil conductivity
        k=ks.*(thetai/thetas).^(2*b+3);
        k1=k(1:m-1);
        k2=k(2:m);
        kavg(2:m)=2*(k1.*k2)./(k1+k2);
        kavg(1)=k1(1);
        %------------------------------------------------------------------  
        % Step 3: Compute distribution of root water uptake and root water
        % potential and leaf water potential based on known soil water 
        % potential distribution and transpiration from previous Steps.
        
        % Root water uptake (RWU)
        [psisB, T_r, Tr_total] = RWU_CW(Al, fe_mm(kk), p, kavg, B, l, Kr, Lr, zs, dzs);
        Tr_layer=T_r; 
        %------------------------------------------------------------------
        % Leaf water potential
        [psiL, Fe] = Leaf_water_potential_CW(fe_mm(kk), Al, psisB); % [Pa]
        %------------------------------------------------------------------        
        % Step 4:
        % Remember save the data of Exfiltration
        % Exfiltration rate (capillary rise): Laio et al., 2009
        % Exfl=max(T_r(N_saturate-1),0);
                
        % Cumulative root water uptake from bottom
        T_r_cumsum=flip(cumsum(flip(T_r)));
      
        % Residual recharge rate [m s-1] accounting for root water uptake
        % in the saturated zone
        Q_re_res=Q_re-T_r_cumsum(N_saturate);
        
        % Compute flux
        gp(2:m)=diff(p)/dz;
        q(2:m)=-kavg(2:m).*gp(2:m);  % flux driven by the states of the neighbourhood cells (theta/psis): [m s^-1]
        
        %----------------- Upper Boundary Condition on Flux
        q(1)=qs_up;
        %----------------- Lower Boundary Condition on Flux----------------------------
        % Zero flux gradient
        q(m+1)=q(m);

        T_r(N_saturate:end)=0; % Make sure T_r is compenseted by the Q_re in the saturated zone (T_r=0)
        if Q_re_res<0  % downward
            if q(N_saturate)<-Q_re_res
                q(1, N_saturate+1:m+1)=-Q_re_res;
                T_r(N_saturate+1:end)=0; 
            else
                q(1, N_saturate:m+1)=-Q_re_res;           
            end
        else
            q(1, N_saturate:m+1)=-Q_re_res;  %upward
        end

        % Compute exfiltration as upward flux
        Exfl=Exfl+min(q(N_saturate-1),0);
        
        % Compute flux gradient
        gq(1:m)=diff(q)./dz;
       
        %-----------------Solve Richards equation with explicit scheme------------
        % (As*dzs)d(theta)/dt=-As*d(q)-(As*dzs)*(B)*Q
        %=>  d(theta)/dt=-d(q)/dzs-(B)*Q
        thetaf(1:m)=thetai(1:m)-(dt)*(gq(1:m))-(dt)*T_r/dz;
           
        % Upward Over flow to neighboring unsaturated cell: ensure conservation
        Q_o=(thetaf(N_saturate-1)-thetas); % Difference between saturation and theta(N_s-1) after updating
        thetaf(N_saturate-1)=min(thetaf(N_saturate-1),thetas); % theta(N_s-1) cannot be larger than saturation
        thetaf(N_saturate-2)=thetaf(N_saturate-2)+max(Q_o,0);  % move excessive water to theta(N_s-2) if happened    
        
             
        %  Make sure saturated zone is saturated (Do not do this)  
        %thetaf(min(find(thetaf==thetas)):end)=thetas;
        
        % Checking water balance
        thetapre=thetai;
        thetai=thetaf;
        V1=V1+(sum(thetaf)-sum(thetapre))*dz;
        V2=V2+(Q_re-T_r_cumsum(1))*dt;
        
    end
    
    % Update averaged psisL for calculating MWUE
    psisL_30avg_oneday= [psisL_30avg_oneday(2:48) psiL/10^6];
    
    % Save data
    z_GT(kk,1)=zw;
    TIME(kk,1)=kk/2;
    QQ(kk,1)=Q_re;
    psi_L(kk,1)=psiL;
    psi_B(kk,1)=psisB*9.8*1000;
    theta_mm(kk,:)=thetaf;
    Tr_layer_mm(kk,:)=Tr_layer;
    Tr_Total_mm(kk,1)=Tr_total;
    EXFL(kk,1)=Exfl;
end
toc
thetaf_sum=sum(thetaf)*dzs;  % Final total water volumn [m3];


DATA.timeseries=table(TIME, z_GT, QQ, psi_L, psi_B, Tr_Total_mm, La_mm, fe_mm, An_mm, gs_mm, EXFL);
DATA.zs=zs;
DATA.thetas=theta_mm;
DATA.Tr_layer=Tr_layer_mm;
save('DATA.mat', 'DATA')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% subplot(1,2,1)
% plot(thetaf,-z,'r-')
% xlabel('\theta_s (m^{3} m^{-3})')
% ylabel('depth (m)')
% 
% subplot(1,2,2)
% plot(p,-z,'r-')
% xlabel('head (m)')
% ylabel('depth (m)')
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% figure
% subplot(1,2,1)
% plot(TIME/24, -z_GT,'k-')
% xlabel('Time (day)')
% ylabel('Groundwater level (m)')
% 
% subplot(1,2,2)
% plot(TIME/24, QQ,'k-')
% xlabel('Time (day)')
% ylabel('Recharge rate (m s^{-1})')
% 
% 
% figure
% subplot(2,1,1)
% hold on
% plot(TIME, psi_L/10^(6),'r-')
% plot(TIME, psi_B/10^(6),'b-')
% ylabel('leaf/root water potential (MPa)','fontweight','bold','fontsize',10)
% subplot(2,1,2)
% plot(TIME, fe_mm,'r-')
% xlabel('Time (day)','fontweight','bold','fontsize',10)
% ylabel('Transpiration rate (mol m^{-2} s^{-1})','fontweight','bold','fontsize',10)
% 
% 
% figure
% subplot(2,1,1)
% plot(TIME, gs_mm,'r-')
% ylabel('Stomatal conductance (mol m^{-2} s^{-1})','fontweight','bold','fontsize',10)
% subplot(2,1,2)
% plot(TIME, An_mm,'r-')
% ylabel('Assimilation rate (mol m^{-2} s^{-1})','fontweight','bold','fontsize',10)
% xlabel('Time (day)','fontweight','bold','fontsize',10)
% 
% figure
% subplot(2,1,1)
% pcolor (TIME/24,-zs,theta_mm')
% colormap(flipud(parula));
% shading ('interp')
% colorbar
% title ('\theta_s (m^{3} m^{-3})','fontweight','bold','fontsize',10)
% xlabel('Time (day)','fontweight','bold','fontsize',10)
% ylabel('z (m)','fontweight','bold','fontsize',10)
% 
% subplot(2,1,2)
% pcolor (TIME/24,-zs,Tr_layer_mm')
% colormap(flipud(parula));
% shading ('interp')
% colorbar
% title ('Root influx (+) or efflux (-) (m s^{-1})','fontweight','bold','fontsize',10)
% xlabel('Time (day)','fontweight','bold','fontsize',10)
% ylabel('z (m)','fontweight','bold','fontsize',10)
% 
% 
% 
% 
% for i=1:length(z_GT)
%     T_GW(i,1)=sum(Tr_layer_mm(i,find(zs>=z_GT(i))));  % T_GW
%     T_US_efflux(i,1)=sum(Tr_layer_mm(i,find(zs<z_GT(i) & Tr_layer_mm(i,:)<0)));  % HR
%     T_US_influx(i,1)=sum(Tr_layer_mm(i,find(zs<z_GT(i) & Tr_layer_mm(i,:)>0)));
% end
% 
% figure
% subplot(3,1,1)
% plot(TIME/24, T_GW)
% ylabel('Root influx_{saturated} (m s^{-1})','fontweight','bold','fontsize',10)
% 
% subplot(3,1,2)
% plot(TIME/24, abs(T_US_efflux))
% ylabel('Root efflux (m s^{-1})','fontweight','bold','fontsize',10)
% 
% subplot(3,1,3)
% plot(TIME/24, T_US_influx)
% xlabel('Time (day)','fontweight','bold','fontsize',10)
% ylabel('Root influx_{unsaturated} (m s^{-1})','fontweight','bold','fontsize',10)
% 
% 


