%% Simulating the excitation/inhibition regulation in Alzheimer's disease with a tripartite synapse model
%% Authors: Shangbin Chen, Jinyu Li 
%% Please contact us for any questions.
%% Email: lijinyu@hust.edu.cn

%% The first glutamate synaptic model: contains the mechanism by which astrocytes are activated by Aβ leading to excitotoxicity


%% Clear Memory & Display Start Time
clear all% Clears memory
Start_Time=clock;% Record start time
fprintf('%s %d %d %d','The Start Time is:',fix(Start_Time(4:6)),' (HH:MM:SS) and')% Display Start Time
tic;% To calculate total time elapsed it will be followed by 'toc' in the 
% last cell.

global V Eca  cin Cout  deltaT  i JN JT JL JR gR gL gT gN F ftf fts Nd Rd Rf Ld Td Lf kVGCC a

%% Time-duration of Simulation
Tmax=5.0e+4;% Time-duration of simulation; unit: ms
deltaT=0.05;% Fixed Time-step; unit: ms
tn=Tmax/deltaT;% Total number of time-points
t=0:deltaT:Tmax;% Time-vector

%% Other Misc. Parameters
tdr=0;% Records previous time at which synaptic vesicle release occured; unit: ms
R=8.314;% Real Gas constant; unit: J/K mole
T = 310;% 273.15+25; unit: K  25℃
F=96487;% Faradays Constant; unit: Coulombs/mole

%% Bouton parameters
rad_neuro=0.3141e-4;% Radius of bouton; unit: cm; Ref: Calculated using vneuro
vneuro=(10^(-3))*(4/3)*pi*((rad_neuro)^3);% Volume of bouton; unit: L; Ref: Koester & Sakmann (2000)
sneuro=4*pi*(rad_neuro)^2;%Surface area of bouton; unit: cm^2; Ref: Koester & Sakmann (2000)

gna=120;%35;% Sodium conductance density; unit: mS/cm^2; Ref: Hodgkin & Huxley (1952)
gk=36;%6;% Pottasium conductance density; unit: mS/cm^2; Ref: Hodgkin & Huxley (1952)
p_ca=2.3e-9;% Conductance of single N-type channel; unit: mS; Ref: Weber et al (2010)
rho_ca=(4*8.0659e+007);% Density of N-type channels; unit: per cm^2; Ref: 
                     % Determined through computer simulations so that the 
                     % average Pr lies between 0.2�0.3 (when astrocyte is not stimulated) 
                     % similar to the experiments of Perea & Araque (2007)
gc=p_ca*rho_ca;% Calcium channel conductance density; unit: mS / cm^2
gl=0.3;% Leak conductance density; unit: mS/cm^2; Ref: Hodgkin & Huxley (1952)
vl=-59.4;% Reversal potential for leak; unit: mV; Ref: Hodgkin & Huxley (1952)
vna=45;% Reversal potential for Sodium ion; unit: mV; Ref: Hodgkin & Huxley (1952)
vk=-82;% Reversal potential for Pottasium ion; unit: mV; Ref: Hodgkin & Huxley (1952)
vca=125;% Reversal potential for Calcium ion determined through Nernst equation 
           % assuming an extracellular calcium concentration of 2 mM which
           % is equal to the extracellular calcium concentration in the
           % experiments of Perea & Araque (2007)
c1=0.185;% Ratio of Volume ER to Volume Neuron cell; unit: dimensionless; Ref: Shuai & Jung (2002)
v1=30*10^(-3);% Maximal calcium flux through IP3R; unit: per ms; Ref: Determined through Simulations
v2=0.2374*10^(-3);% Leak of calcium from ER to cytosol; unit: per ms; Ref: Determined through Simulations
v3=90;% Maximal SERCA pump rate; unit: nM/ms; Ref: Determined through Simulations
k3=0.1e+3;% Michaelis-Menten (MM) constant for SERCA; unit: nM; Ref: Erler et al (2004)
Ip=0.4;% Maximum PMCa current; unit: uA per cm2; Ref: Determined through computer simulations
k_pump=100;% MM constant for PMCa; unit: nM; Ref: Erler et al (2004)
v_leak=0.001022664392140;% Maximum leak of Ca2+ through plasma membrane; unit: per ms; Ref: Determined through simulations
%IP3R-kinetics 
d1=0.13*10^(3);% IP3 dissociation constant; unit: nM; Ref: Shuai & Jung (2002)
d2=1.049*10^(3);% Inhibitory Ca2+ dissociation constant; unit: nM; Ref: Shuai & Jung (2002)
d3=0.9434*10^(3);% IP3 dissociation constant; unit: nM; Ref: Shuai & Jung (2002)
d5=0.08234*10^(3);% Activation Ca2+ dissociation constant; unit: nM; Ref: Shuai & Jung (2002)
a2=0.2*10^(-6);% Inhibitory Ca2+ binding constant; unit: per nM ms; Ref: Shuai & Jung (2002)
% Calcium current parameters
v_half=-17;% Half-activation voltage; unit: mV; Ref: Ishikawa et al (2005)
kca=8.4;% Slope factor; unit: mV; Ref: Ishikawa et al (2005)
taumc=10;% Time-constant; unit: ms; Ref: Ishikawa et al (2005)
%IP3 generation parameters
v_glu=0.062;% Maximal IP3 production rate from mGluRs; unit: nM per ms; Ref: Nadkarni & Jung (2008)
k_glu=0.78e-3;% Glutamate concentration at which v_glu is halved; unit: mM; Ref: Nadkarni & Jung (2008)
tau_ip3=0.14e-3;% IP3 degradation constant; unit: per ms; Ref: Nadkarni & Jung (2008)
np=0.3;
% Vesicle Regulatory time constants; unit: ms; Ref: Tsodyks & Markram (1998)
tau_rec=800;% Vesicle recovery time constant; unit: ms; Ref: Tsodyks & Markram (1997)
tau_inact=3;% Vesicle inactivation time constant; unit: ms; Ref: Tsodyks & Markram (1997)

vr=20;% Resting membrane potential of bouton; unit: mV
%gv=1;%60;% Glutamate concentration in single vesicle; unit: mM; Ref: Montana et al. (2006)

% Synaptic vesicle parameter
Omega_c=40e-3;% Neurotransmitter clearance rate s-1
rho_c=0.005;%synaptic vesicle-to-extracellular space volume ratio
Y_T=100;%500;%.*mmole
omega_f=2e-3;% ms-1
omega_d=2e-3;% ms-1
U0=0.15;

% Influence parameters of glial transmitter on the probability of presynaptic glutamate release
ksi=0.5;%  the parameter lumps,in a phenomenological way, the information on the effect of gliotransmission on synaptic release.
O_p=1;%  the rise rate of the effect of gliotransmission on synaptic glutamate release
tau_p=120;% the decay time of the effect of gliotransmission on synaptic glutamate release
%% Spine Parameters
R_in=0.7985e+8;% Input resistance of dendrite spine; unit: k-ohm; Ref: Calculated
tau_mem=50;% Post-synaptic membrane time constant; unit: ms; Ref: Tsodyks & Markram (1997)
rad_spine=0.6e-4;% Radius of spine-head; unit: cm; Ref: Dumitriu et al (2010)
vspine=(10^(-3))*(4/3)*pi*((rad_spine)^3);% Volume of spine-head assuming it to be spherical in shape; unit: L; Ref: Calculated
sspine=4*pi*(rad_spine)^2;% Surface area of spine-head assuming it to be spherical in shape; unit: cm^2; Ref: Calculated
v_post=-70;% Resting membrane potential of spine membrane; unit: mV
ks=100e-3;% Calcium extrusion rate by PMCa; unit: per ms; Ref: Keller et al (2008)

% HH model
gna_post=45;%35;% Sodium conductance density; unit: mS/cm^2; Ref: 2018
gk_post=18;%6;% Pottasium conductance density; unit: mS/cm^2; Ref: 2018
gl_post=0.05;% Leak conductance density; unit: mS/cm^2; Ref: 2018
vna_post=55;% Reversal potential for Sodium ion; unit: mV; Ref: 2018
vk_post=-90;% Reversal potential for Pottasium ion; unit: mV; Ref: 2018
vl_post=-70;% Reversal potential for leak; unit: mV; Ref: 2018

% AMPA Receptor Parameters 
alpha_ampa=1.1;% Binding of transmitter to receptor; unit: mM per ms
beta_ampa=0.19;%0.67;% Dissociation of transmitter from receptor; unit: per ms
g_ampa=0.035;%0.0145;%0.35e-6; Conductance of AMPAR; unit: mS

% NMDA Receptor Parameters 
alpha_nmda=0.072;%2.2;% Binding of transmitter to NMDAR; unit: mM/ms
beta_nmda=0.066;%0.67;% Dissociation of transmitter from NMDAR; unit: ms-1
g_nmda=0.035;%0.026;%0.2e-6;% Conductance of NMDAR; unit: mS
g_nmdapre=0.1;
bt=200e+3;% Endogenous buffer concentration; unit: nM;
K_bt=10e+3;% Ca2+ affinity of Endogenous buffer; unit: nM

%% Astrocyte Parameters
% All values are from De Pitta et al (2009) except 'aNa.'
aNa=20;% Number of IP3 channels in cluster; unit: dimensionless; Ref: Nadkarni & Jung (2007)
ac0=2e+3;% Total cell free Ca2+ concentration; unit: nM
ac1=0.185;% Ratio of ER volume to cytosol volume; unit: dimensionless
av1=6e-3;% Maximal IP3R flux; unit: per ms
av2=0.11e-3;% Maximal rate of Ca2+ leak from ER; unit: per ms  from Bronac Flanagan
av3=2.2;%0.9;% Maximal rate of SERCA uptake; unit: nM per ms  
ak3=0.05e+3;%0.1e+3;% SERCA Ca2+ affinity; unit: nM
ad1=0.13e+3;% IP3 dissociation constant; unit: nM
ad2=1.049e+3;% Ca2+ inactivation dissociation constant; unit: nM
ad3=0.9434e+3;% IP3 dissociation constant; unit: nM
ad5=0.08234e+3;% Ca2+ activation dissociation constant; unit: nM
aa2=0.2e-6;% IP3R binding rate for Ca2+ Inhibition; unit: per nM ms

Vast=3.49*1e-13; % Volume of an astrocyte, unit: liter

%RyR  Jin pm  para          
k0=0.005e-3;%0.013e-3;% Zero calcium concentration level leak from RyRs unit: ms-1 frmo Liu
k2=0.18e-12;%0.18e-12; % Maximal rate of the RyRs unit: ms-1 frmo Liu
kd=0.13e+3;%0.13e+3; % RyR sensitivity for the CICR unit: nM frmo Liu
k1=0.95e-3; %0.5e-3; % Rate constant of calcium extrusion unit: ms-1 frmo Liu

%VGCC para
Cout=1500; % Ca2+ concentration in ECS, unit: uM
gT=0.0600; % Steady conductance of T type channel, unit: pS
gL=3.5000; % Steady conductance of L type channel, unit: pS
gN=0.3900; % Steady conductance of N type channel, unit: pS
gR=0.2225; % Steady conductance of R type channel, unit: pS
Nl=1;        %%% open or block the L type channel,  1=open, 0=close 
Nn=1;        %%% open or block the N type channel,  1=open, 0=close
Nr=1;        %%% open or block the R type channel,  1=open, 0=close
Nt=1;        %%% open or block the T type channel,  1=open, 0=close
 
%TRPA1
delta_H = -28;%kJ/mol
delta_S = -107;%J/(mol k)
gTRPA1 = 2;%pS
N = 1;%from 2017 Marie Mulier
Ltrpa = 4.4e-4;%from 2011 Marcelo Salazar
Ctrpa = 1e4;%from 2011 Marcelo Salazar
Ftrpa = 5e3;%0.11;%from 2011 Marcelo Salazar
K_D = 20;%0.13;%uM
Ligand = 1;%uM
Qtrpa = Ligand/K_D;
Ktrpa = exp(-(delta_H-T*delta_S)/(R*T));
i_to_j = 1e4/(3.49)/(2*F);
V_rev = 0;%from "Allyl isothiocyanate (AITC) activates nonselective cation currents in human cardiac fibroblasts: possible involvement of TRPA1"
          %or "Transient Receptor Potential Channel A1 Is Directly Gated by Calcium Ions*"

% IP3 regulation parameters
av_plcd=0.05;% Maximal rate of IP3 production by PLC?; unit: nM per ms 
ak_plcd=1.5e+3;% Inhibition constant of PLC? activity; unit: nM
aK_plcd=0.1e+3;% Ca2+ affinity of PLC?; unit: nM
ar5p=0.05e-3;% Maximal rate of degradation by IP-5P; unit: per ms
av_plcb=0.05;% Maximal rate of IP3 production by PLC?; unit: nM per ms
aK_R=1.3e-3;% Glutamate affinity of the receptor; unit: mM
aK_P=10e-3;% Ca2+/PKC-dependent inhibition factor; unit: mM
aK_pi=0.6e+3;% Ca2+ affinity of PKC; unit: nM
av_3K=2;% Maximal rate of degradation by IP3-3K; unit: nM per ms
aK_D=0.7e+3;% Ca2+ affinity of IP3-3K; unit: nM
aK3=1e+3;% IP3 affinity of IP3-3K; unit: nM

% Glutamate regulatory parameters
degG=0.02;%0.12;% Extracellular glutamate degradation rate; unit: per ms; Ref: Destexhe et al (1998)

% Gliotransmitter release parameters
C_Theta=196.9; %Ca^2+ threshold for exocytosis unit:*umole
Omega_A=0.6e-3; %Gliotransmitter recycling rate unit:per ms
U_A=0.6;     %Gliotransmitter release probability
G_T= 500;%1500;200;%Total vesicular gliotransmitter unit:*mmole
rho_e=6.5e-4; %Ratio of astrocytic vesicle volume/ESS volume
Omega_e=5e-3;%5e-3; %Gliotransmitter clearance rate (think about distributed release)unit:per ms

%Influence of Abeta on calcium channel
kVGCC=10; % The strength of the influence of Abeta on VGCCs
%kin=1;  The strength of the influence of Abeta on abeta channels
kRyR=0.2; % The strength of the influence of Abeta on RyRs
kPLCb=0.05; % The strength of the influence of Abeta on glutamate-dependent IP3 production
kPLCd=0.5; % The strength of the influence of Abeta on Ca2+-dependent IP3 production

%% FIGURE
spk=[]; % Frequency series 
freqs=[]; % Frequency series 
mps=[]; % Abeta levels series 
count=0; % For counting
data=[];
figure
set(gcf,'unit','centimeters','position',[10 10 23 17])
% for aaa=0:1:100
%     a=aaa/100;
for a=[0.2 0.8]
    Ligand = 1+a*5;
    Qtrpa = Ligand/K_D;

%% Variable Initialization
temp=zeros(1,tn+1);% Temporary array
Iapp=zeros(1,tn);% Applied current density
% Pre-synaptic Bouton Variables
G_syn=temp;% Synaptic glutamate concentration
u_syn=temp;
x_syn=temp;
c=temp;% Calcium concentration
cer=temp;% ER Calcium concentration
v=temp;% Membrane potential
m=temp;% Sodium channel activation 
h=temp;% Sodium channel inactivation
n=temp;% Pottassium channel activation
p=temp;% IP3 concentration
q=temp;% IP3 gating variable
mc=temp;% VGCC gating variable
R_syn=temp;% Fraction of releasable vesicles
E_syn=temp;% Fraction of effective vesicles in synaptic cleft
I_syn=temp;% Fraction of inactivated vesicles 
RRP=zeros(1,tn);% Vesicles ready to be released (spontaneous or evoked) 
newRRP=zeros(1,tn);%加了胶质递质对突触释放频率影响后的囊泡

am_nmdapre=temp;
aI_nmdapre=temp;
gama=temp;
mpre_gabaa=temp;% GABAAR gating variable
Ipre_gabaa=temp;% Current through GABAAR

% Astrocyte Variables
ca=temp;% Calcium concentration
auer=temp;% ER Calcium concentration
ax=temp;% IP3R gating variable
ap=temp;% IP3 concentration
aO1=temp;% S1 with calcium bound 
aO2=temp;% S2 with calcium bound
aO3=temp;% S3 with calcium bound
aE_syn=temp;% Fraction of effective SLMV in extra-synaptic cleft 
aI_syn=temp;% Fraction of inactivated SLMV
aR_syn=temp;% Fraction of releasable SLMV
aG_syn=temp;% Glutamate concentration in extra-synaptic cleft

%Astrocyte VGCC
cin=temp;% Calcium concentration weimi
Vin=temp;        %%% Total Ca flux via VGCCs, unit: nM
Ld=temp;         %%% Ld is activation variable in L type channel
Lf=temp;         %%% Lf is inactivation variable in L type channel
Nd=temp;         %%% Nd is activation variable in N type channel
Td=temp;         %%% Td is activation variable in N type channel
Rd=temp;         %%% Rd is activation variable in R type channel
Rf=temp;         %%% Rf is inactivation variable in R type channel
ftf=temp;        %%% ftf is part-variable of inactivation in T type channel
fts=temp;        %%% fts is part-variable of inactivation in T type channel
JN=temp;         %%% JN is the Ca current via N type channel 
JL=temp;         %%% JL is the Ca current via N type channel 
JR=temp;         %%% JR is the Ca current via N type channel 
JT=temp;         %%% JT is the Ca current via N type channel 
V(1:tn)=-75;      %% membrane potential (mV), ranged from -75 to -60 mV 

%Astrocyte GABA
Va(1:tn)=-80.00e-3;
GABA_ast=temp;%unit:M 
GABA_syn=temp;%unit:M
Na_ast=temp;%unit:M
Na_syn=temp; %unit:M
K_ast=temp;%unit:M
K_syn=temp;%unit:M
Glu_ast=temp;%unit:M
% Glu_syn=temp;%unit:M
Ca_syn=temp;%unit:M
x_A=temp;

%Astrocyte TRPA1
Po = temp;%from 2011 Marcelo Salazar
aiTRPA1 = temp;%from 2011 Marcelo Salazar
ajTRPA1 = temp; 

% Post-synaptic spine variables
V_post=temp;% Membrane potential
c_post=temp;
m_ampa=temp;% AMPAR gating variable
I_post=temp;% Current through AMPAR
m_nmda=temp;% NMDAR gating variable
I_nmda=temp;% Current through NMDAR
m_gabaa=temp;% GABAAR gating variable
I_gabaa=temp;% Current through GABAAR
m_info=temp;
b_post=temp;
m_post=temp;% Sodium channel activation 
h_post=temp;% Sodium channel inactivation
n_post=temp;

am_ampa=temp;
aI_post=temp;
am_nmda=temp;
aI_nmda=temp;
I_sic=temp;
rr=temp;

%% Initial conditions
v(1)=-70;% Resting membrane potential of bouton; unit: mV
m(1)=0.1;% Gating variable for sodium channel (activation)
h(1)=0.6;% Gating variable for sodium channel (inactivation)
n(1)=0.3;% Gating variable for potassium channel (activation)
p(1)=160;% Resting [IP3]; unit: nM
c(1)=100;% Resting [Ca^{2+}]; unit: nM
cer(1)=(1.30e+6);% Resting ER [Ca^{2+}]; unit: nM
q(1)=0.22;% Gating variable for IP3R
mc(1)=0;% Gating variable for Calcium channel
mu=[1 0 0 0 0 0];% Inital vector for vesicle '1' & '2.' 
R_syn(1)=1;% Relasable fraction of vesicles
E_syn(1)=0;% Effective fraction of vesicles in synaptic cleft
I_syn(1)=0;% Inactivated fraction of vesicles
G_syn(1)=0.002;% Glutamate concentration in cleft; unit: mM
u_syn(1)=0.8;
x_syn(1)=0.2;
% Peri-synaptic Astrocyte variables
ca(1)=100;% Resting calcium concentration; unit: nM
auer(1)=7000;
ax(1)=0.5;% Gating variable for IP3R 
aO1(1)=0.48;% 1st calcium binding site of SLMV
aO2(1)=0.2;% 2nd calcium binding site of SLMV
aO3(1)=5e-4;% 3rd calcium binding site of SLMV
aE_syn(1)=0;% Fraction of effective SLMV in extra-synaptic cleft
aI_syn(1)=0;% Fraction of inactivated SLMV
aR_syn(1)=1;% Fraction of releasable SLMV
ap(1)=160;% Resting [IP3]; unit: nM
aG_syn(1)=1e-3;% Basal glutamate in the extra-synaptic cleft; unit: mM
x_A(1)=1;
% Post-synaptic spine variables
V_post(1)=-70;% Resting post-membrane potential; unit: mV
m_ampa(1)=0;% Gating variable of AMPAR
I_post(1)=0;% Current through AMPAR; unit: uA
c_post(1)=100;
m_nmda(1)=0;
I_nmda(1)=0;
m_gabaa(1)=0;
I_gabaa(1)=0;
m_info(1)=0;
b_post(1)=0;
m_post(1)=0.1;% Gating variable for sodium channel (activation)
h_post(1)=0.6;% Gating variable for sodium channel (inactivation)
n_post(1)=0;%0.3;

am_ampa(1)=0;
aI_post(1)=0;
am_nmda(1)=0;
aI_nmda(1)=0;

% pre-synaptic botton var
mpre_gabaa(1)=0;
Ipre_gabaa(1)=0;
am_nmdapre(1)=0;
aI_nmdapre(1)=0;
gama(1)=0;

Vin(1)=0.00;     %%% initial Total Ca flux via VGCCs, unit: uM  


for i=1:tn
    %% Applied Current Frequency
% Regular Spike: Here 'mod' is the standard modulus function. 't' is 
% in 'ms.' At present the spike frequency is 5 Hz with duration of 4 ms.    
    if   mod(t(i),200)>=0 && mod(t(i),200)<=2 %t(i)>=20000 &&
      Iapp(i)=5;% Applied current density
    else
      Iapp(i)=0;
    end
    
    %% Bouton Processes
    % Gating Variables
    an=0.01*((-v(i)-60)/(exp((-v(i)-60)/10)-1));bn=0.125*exp((-v(i)-70)/80);% Opening closing parameters for activation of potassium channel
    am=0.1*((-v(i)-45)/(exp((-v(i)-45)/10)-1));bm=4*exp((-v(i)-70)/18);% Opening closing parameters for activation of sodium channel opening
    ah=0.07*exp((-v(i)-70)/20);bh=1/(exp((-v(i)-40)/10)+1);% Opening closing parameters for inactivation of sodium channel opening
    aq=a2*d2*(p(i)+d1)/(p(i)+d3);bq=a2*c(i);% Opening and closing parameters for IP3R gating
    mcinf=1/(1+exp((v_half-v(i))/kca));% VGCC current activation function
    % Ionic Currents
    ina=gna*(v(i)-vna);% Sodium current; unit: uA per cm^2
    ik=gk*(v(i)-vk);% Potassium current; unit: uA per cm^2
    il=gl*(v(i)-vl);% Leak current; unit: uA per cm^2
    ica=gc*(mc(i)^2)*(v(i)-vca);% VGCC current; unit: uA per cm^2
    I_pump=Ip*(c(i)^2/(c(i)^2+k_pump^2));% PMCa current; unit: uA per cm^2
    I_leak=v_leak*(v(i)-vca);% Calcium leak current; unit: uA per cm^2

    % IP3R Kinetics
    minf=p(i)/(p(i)+d1);ninf=c(i)/(c(i)+d5);% Steaty-state functions governing IP3 & Ca2+ based closing of IP3R
    jchan=c1*v1*(minf^3)*(ninf^3)*(q(i)^3)*(c(i)-cer(i));% Calcium flux through IP3R
    jpump=v3*(c(i)^2)/(k3^2+c(i)^2);% Serca pump flux
    jleak=c1*v2*(c(i)-cer(i));% Calcium leak through ER into cytosol
    p_glu=v_glu*((aG_syn(i)^np)/(k_glu^np + aG_syn(i)^np));% IP3 production due to extra-synaptic glutamate
    
    %% Astrocyte Processes
    %IP3 production & degradation terms. All are from De Pitta et al (2009)
    aplcb_ca=1+(aK_P/aK_R)*(ca(i)/(ca(i)+aK_pi));% Calcium-dependent inhibition of 'ap_glu'
    ap_glu=(av_plcb+kPLCb*a) *(G_syn(i)^0.7)/((G_syn(i)^0.7)+(aK_R*aplcb_ca)^0.7);% Agonist-dependent IP3 production
    ap_plcd=(av_plcd+kPLCd*a) *(1/(1+(ap(i)/ak_plcd)))*(ca(i)^2)/((ca(i)^2)+aK_plcd^2);% Agonist-independent IP3 production这里和师兄的公式不太一样
    ap_mapk=av_3K*((ca(i)^4)/((ca(i)^4)+aK_D^4))*(ap(i)/(ap(i)+aK3));% IP3 degradation due to IP3-3K
    ap_deg=ar5p*(ap(i));% IP3 degradation due to IP-5P
    
    % IP3R gating variables
    aaq=aa2*ad2*(ap(i)+ad1)/(ap(i)+ad3);abq=aa2*ca(i);% gating variables for IP3R
    aminf=ap(i)/(ap(i)+ad1);aninf=ca(i)/(ca(i)+ad5);% Steaty-state functions governing IP3 & Ca2+ based closing of IP3R    
%     auer=(ac0-ca(i))/ac1;% ER Calcium concentration

    % Using Box-Muller Algorithm (Fox, 1997) to simulate noise-term in 
    % equation (12) of the manuscript
    au1=rand(1,1);au2=rand(1,1);% Two uniformaly distributed random variables
    aa=((aaq*(1-ax(i))-abq*ax(i)))/aNa;% Calculating co-variance at each time-step 
    dW=sqrt(-(2*deltaT*(aa)*log(au1)))*cos(2*pi*au2);% Generating independent Gaussian random number
    % Because Gaussian random numbers can have arbitrarily large positive 
    % or negative values, it is possible that in the course of the 
    % simulation a ax value will fall outside the range of [0,1].Thus a
    % restriction has been imposed to chose those values of dW for which ax 
    % lies in the range [0,1] and discard all other values.
    if ax(i)>=0 && ax(i)<=1 
        ax(i+1)=deltaT*(aaq*(1-ax(i))-abq*ax(i))+ax(i)+dW;
    else
        ax(i+1)=deltaT*(aaq*(1-ax(i))-abq*ax(i))+ax(i);
    end
   
    ajchan=av1*(aminf^3)*(aninf^3)*((abs(ax(i)))^3)*(ca(i)-auer(i));% Calcium flux through IP3R
    ajpump=av3*(ca(i)^2)/(ak3^2+ca(i)^2);% SERCA pump flux
    ajleak=av2*(ca(i)-auer(i));% Calcium leak from ER into cytosol  
    ajRyR=(k0+((k2*ca(i)^3)/( (kd+kRyR*a)^3 +ca(i)^3)))*(ca(i)-auer(i));  % from Liu
    ajpm=k1*ca(i); % from Liu
    %ajin=v5+kin*a^m; % from Liu

    Po(i) = 1/(1+(1+Ktrpa+Qtrpa+Ktrpa*Qtrpa)/(Ltrpa*(1+Ktrpa*Ctrpa+Qtrpa*Ftrpa+Ktrpa*Ctrpa*Qtrpa*Ftrpa)));%from 2011 Marcelo Salazar
    aiTRPA1(i) = N*gTRPA1*Po(i)*(V(i)-V_rev);%from 2011 Marcelo Salazar
    ajTRPA1(i) = -i_to_j*0.17*aiTRPA1(i);   
    
    cin(i)=ca(i)*0.001;
    Eca(i)=59.5/2*log10(Cout/cin(i)); %%% Eca is Nernst Potential of Ca, 59.5=R*T/F*1000*ln10, unit as mV
    NRTC();           %% to get the Ca current via R type channel, JR  
    NLTC();           %% to get the Ca current via L type channel, JR 
    NNTC();           %% to get the Ca current via N type channel, JR 
    NTCC();           %% to get the Ca current via T type channel, JR 
    in(i)=(Nn*JN(i)+Nl*JL(i)+ Nr*JR(i)+Nt*JT(i));    %%% Total Ca current via VGCCs %%%%%%%%%%%%%%%%%%% 
    Vin(i)=-in(i)*1e4/(3.49)/(2*F);                 %%% Converting the current to the flux with unit as nM/ms  
 

    %% Spine Processes
    % Gating Variables
    b_info=1/(1+exp((V_post(i)+45)/7));
    tau_b=0.1+0.75/(1+exp( -(V_post(i)+40.5) /(-6) ));
    n_info=1/(1+exp(-(V_post(i)+35)/10));
    tau_n=0.1+0.5/(1+exp((V_post(i)+27)/15));

%     an_post=0.01*((-V_post(x)-60)/(exp((-V_post(x)-60)/10)-1));bn_post=0.125*exp((-V_post(x)-70)/80);% Opening closing parameters for activation of potassium channel
%     am_post=0.1*((-V_post(x)-45)/(exp((-V_post(x)-45)/10)-1));bm_post=4*exp((-V_post(x)-70)/18);% Opening closing parameters for activation of sodium channel opening
%     ah_post=0.07*exp((-V_post(x)-70)/20);bh_post=1/(exp((-V_post(x)-40)/10)+1);% Opening closing parameters for inactivation of sodium channel opening
    % Ionic Currents
    ina_post=gna_post*(V_post(i)-vna_post);% Sodium current; unit: uA per cm^2
    ik_post=gk_post*(V_post(i)-vk_post);% Potassium current; unit: uA per cm^2
    il_post=gl_post*(V_post(i)-vl_post);% Leak current; unit: uA per cm^2

    Mg=1/(1+exp(-0.062*V_post(i))*2/3.57);
    I_post(i)=g_ampa*m_ampa(i)*(V_post(i));% AMPAR current; unit: uA
    I_nmda(i)=g_nmda*Mg*m_nmda(i)*V_post(i);
    aI_post(i)=g_ampa*am_ampa(i)*V_post(i);% AMPAR current; unit: uA
    aI_nmda(i)=g_nmda*Mg*am_nmda(i)*V_post(i);
    I_sic(i)=aI_post(i)+aI_nmda(i);
    %% Bouton DF
    p(i+1)=p(i)+deltaT*(((p(1)-p(i))*(tau_ip3)));% Intracellular [IP3]  +p_glu
    m(i+1)=m(i)+deltaT*(am*(1-m(i))-bm*m(i));% Sodium channel activation
    h(i+1)=h(i)+deltaT*(ah*(1-h(i))-bh*h(i));% Sodium channel inactivation
    n(i+1)=n(i)+deltaT*(an*(1-n(i))-bn*n(i));% Potassium channel activation
    q(i+1)=q(i)+deltaT*(aq*(1-q(i))-bq*q(i));% IP3R gating variable
    mc(i+1)=mc(i)+deltaT*((mcinf-mc(i))/taumc);% VGCC gating variable

    v(i+1)=v(i)+deltaT*(Iapp(i)-((m(i)^3)*h(i)*ina+(n(i)^4)*ik+il));% Membrane potential displacemnt from rest   +aI_nmdapre(i)
    c(i+1)=c(i)+deltaT*(((-ica-I_pump-I_leak)*sneuro)/(2*F*vneuro)-jpump-jleak-jchan);% Intracellular calcium 
    cer(i+1)=cer(i)+deltaT*((jchan+jleak+jpump)/c1);% ER calcium concentration
    
    gama(i+1)= gama(i)+deltaT*((1-gama(i))*O_p*aG_syn(i)-(gama(i)/tau_p));
    u_syn(i+1)=u_syn(i)+deltaT*((-u_syn(i))*omega_f);
    x_syn(i+1)=x_syn(i)+deltaT*( (1-x_syn(i))*omega_d);
    G_syn(i+1)=G_syn(i)+deltaT*( -degG*(G_syn(i)));%(I_Gluast(i)*S_A/(Z_Glu*F*Vol_S))

    % check if spiking
    if  v(i+1)>vr && v(i)<vr    % The following formulation makes sure that a "spike" is only triggered at the first threshold crossing
            U=U0+(ksi-U0)*gama(i);%
            u_syn(i+1)=u_syn(i)+ U* (1-u_syn(i)) ;
            r_syn=u_syn(i+1) * x_syn(i);
            x_syn(i+1)=x_syn(i) - r_syn ;
            G_syn(i+1)=G_syn(i)+ rho_c * Y_T * r_syn;
            flag=1;
    end
    
    %% Spine DF
    m_info(i+1)=1/(1+exp(-(V_post(i)+30)/9.5));
    b_post(i+1)=b_post(i)+deltaT*((b_info-b_post(i))/tau_b);
    n_post(i+1)=n_post(i)+deltaT*((n_info-n_post(i))/tau_n);

    V_post(i+1)=V_post(i)+deltaT*( -( (m_info(i)^3)*b_post(i)*ina_post+(n_post(i)^4)*ik_post+il_post+I_post(i)+I_nmda(i)+I_sic(i) ) );%+I_gabaa(i)

    theta=bt*K_bt/((K_bt+c_post(i))^2);
    m_ampa(i+1)=m_ampa(i)+deltaT*(alpha_ampa*G_syn(i)*(1-m_ampa(i))-beta_ampa*m_ampa(i));% AMPAR gating variable
    m_nmda(i+1)=m_nmda(i)+deltaT*(alpha_nmda*(G_syn(i))*(1-m_nmda(i))-beta_nmda*m_nmda(i));
    am_ampa(i+1)=am_ampa(i)+deltaT*(alpha_ampa*aG_syn(i)*(1-am_ampa(i))-beta_ampa*am_ampa(i));% AMPAR gating variable
    am_nmda(i+1)=am_nmda(i)+deltaT*(alpha_nmda*aG_syn(i)*(1-am_nmda(i))-beta_nmda*am_nmda(i));
    
    c_post(i+1)=c_post(i)+deltaT*((-(0.012*I_post(i)+0.12*I_nmda(i))/(2*F*vspine)-ks*(c_post(i)-c_post(1)))/(1+theta));

    %% Astrocyte DF
    ca(i+1)=ca(i)+(-ajchan-ajpump-ajleak-ajRyR-ajpm+Vin(i)+ajTRPA1(i))*deltaT;% Dynamic calcium concentration  
    
    ap(i+1)=ap(i)+deltaT*(ap_plcd-ap_mapk-ap_deg);% IP3 concentration  +ap_glu
    auer(i+1)=auer(i)+deltaT*((ajpump+ajchan+ajleak+ajRyR)/ac1);  % ER Calcium concentration  

    x_A(i+1)=x_A(i)+deltaT*(Omega_A * (1 - x_A(i)));%Fraction of gliotransmitter resources available for release
    aG_syn(i+1)=aG_syn(i)+deltaT*(-Omega_e * aG_syn(i));%gliotransmitter concentration in the extracellular space

    % check if spiking
    if ca(i+1)>C_Theta && ca(i)<C_Theta     % The following formulation makes sure that a "spike" is only triggered at the first threshold crossing
            aG_syn(i+1) = aG_syn(i) + rho_e * G_T * U_A * x_A(i);%The gliotransmitter release happens when the threshold is crossed, in Brian terms it can therefore be considered a "reset"
            x_A(i+1) = x_A(i) - U_A *  x_A(i);
    end
    
end

%Store key result parameters
data = [data;ca;G_syn;aG_syn;V_post];

%% Plot a contrast diagram of the results at different Aβ concentrations
subplot(4,2,1+count),plot(ca,'k','Linewidth',2)
set(gca,'xLim',[0,Tmax/deltaT],'FontName','Arial','FontSize',12,'LineWidth',1.5,'Fontweight','bold','Box','off')
xticks(linspace(0,length(ca)-1,3))
xticklabels('')
ylabel('[Ca^{2+}](nM)','FontName','Arial','FontSize',14)
title('Astrocyte Ca^{2+}','FontName','Arial','Fontsize',14,'Fontweight','bold')

hold on
subplot(4,2,3+count),plot(G_syn,'k','Linewidth',2)
set(gca,'xLim',[0,Tmax/deltaT],'FontName','Arial','FontSize',12,'LineWidth',1,'Fontweight','bold','Box','off')
xticks(linspace(0,length(G_syn)-1,3))
xticklabels('')
ylabel('g(mM)','FontName','Arial','FontSize',14)
title('Synaptic glutamate','FontName','Arial','Fontsize',14,'Fontweight','bold')
ylim([0 0.3])

hold on
subplot(4,2,5+count),plot(aG_syn,'k','Linewidth',2)
set(gca,'xLim',[0,Tmax/deltaT],'FontName','Arial','FontSize',12,'LineWidth',2,'Fontweight','bold','Box','off')
xticks(linspace(0,length(aG_syn)-1,3))
xticklabels('')
ylabel('g(mM)','FontName','Arial','FontSize',14)
title('Extra-synaptic glutamate','FontName','Arial','Fontsize',14,'Fontweight','bold')

hold on
subplot(4,2,7+count),plot(V_post,'k','Linewidth',2)
set(gca,'xLim',[0,Tmax/deltaT],'FontName','Arial','FontSize',12,'LineWidth',2,'Fontweight','bold','Box','off')
xticks(linspace(0,length(V_post)-1,3))
xticklabels('')
xticklabels(linspace(0,length(V_post)-1,3)/20000)
ylabel('V(mV)','FontName','Arial','FontSize',14)
xlabel('Time(sec)','FontName','Arial','FontSize',14)
title('Postsynaptic membrane potential','FontName','Arial','Fontsize',14,'Fontweight','bold')
hold on
count=count+1;
end


toc;% To display total time elapsed


%%%%%%%%%%%% Different  VGCCs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% L-type %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NLTC()   %%%%% L type channel, JL is the Ca current.
global V deltaT cin i Eca gL Ld JL Lf kVGCC a
if i==1
    Ld(1)=1/(1+exp(-(V(i)+50)/3));  
    return;
else
   Ld(i)=Ld(i-1);
end;
    und(i)=1/(1+exp(-(V(i)+50)/3)); 
    taod(i)=18*exp(-((V(i)+45)/20)^2)+1.5;
    Lf(i)=0.00045/(0.00045+cin(i)/1000);
    Ld(i)=(und(i)-Ld(i))*deltaT/taod(i)+Ld(i-1);  
    dd=Ld(i);
    Ld(i)=(und(i)-Ld(i))*deltaT/taod(i)+Ld(i-1);
while((sum(dd-Ld(i))>0.01)||(sum(dd-Ld(i))<-0.01))
    dd=Ld(i);
    und(i)=1/(1+exp(-(V(i)+50)/3));   
    taod(i)=18*exp(-((V(i)+45)/20)^2)+1.5;
    Lf(i)=0.00045/(0.00045+cin(i)/1000);
    Ld(i)=(und(i)-Ld(i))*deltaT/taod(i)+Ld(i-1);
 end;
    JL(i)=(gL+kVGCC*a)*Ld(i).*Lf(i).*(V(i)-Eca(i));  % JL is the Ca current, unit: pA
end

%%% N-type %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NNTC()   %%%%% N type channel, JN is the Ca current
global V gN deltaT und cin  i   Eca JN Nd

if i==1
    Nd(1)=1/(1+exp(-(V(i)+45)/7));
    return;
else
    Nd(i)=Nd(i-1);
end;
 und(i)=1/(1+exp(-(V(i)+45)/7));
   
    taod(i)=18*exp(-((V(i)+70)/25)^2)+0.3;
    Nf(i)=0.0001/(0.0001+cin(i)/1000);
    Nd(i)=(und(i)-Nd(i))*deltaT/taod(i)+Nd(i-1);  
    dd=Nd(i);
    Nd(i)=(und(i)-Nd(i))*deltaT/taod(i)+Nd(i-1);  
while((sum(dd-Nd(i))>0.01)||(sum(dd-Nd(i))<-0.01))
        dd=Nd(i);
    und(i)=1/(1+exp(-(V(i)+45)/7));   
    taod(i)=18*exp(-((V(i)+70)/25)^2)+0.3;
    Nf(i)=0.0001/(0.0001+cin(i)/1000);
    Nd(i)=(und(i)-Nd(i))*deltaT/taod(i)+Nd(i-1);    
end;
    JN(i)=gN*Nd(i).*Nf(i).*(V(i)-Eca(i));  %%JN is the Ca current via N type channel
end

%%% N-type %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function NRTC() %%%%% R type channel, JR is the Ca current.
global gR V deltaT und  Rd Rf i  dd  Eca JR ff unf 
if i==1;
  Rd(1)=1/(1+exp(-(V(i)+10)/10));
  Rf(1)=1/(1+exp((V(i)+48)/5));
  return;
else
  Rd(i)=Rd(i-1);
  Rf(i)=Rf(i-1);
end
    und(i)=1/(1+exp(-(V(i)+10)/10));
    unf(i)=1/(1+exp((V(i)+48)/5));
    taod(i)=0.1*exp(-((V(i)+62)/13)^2)+0.05;
    taof(i)=0.5*exp(-((V(i)+55.6)/18)^2)+0.5;
    Rd(i)=(und(i)-Rd(i))*deltaT/taod(i)+Rd(i-1);
    dd=Rd(i);
    Rf(i)=(unf(i)-Rf(i))*deltaT/taof(i)+Rf(i-1);
    ff= Rf(i);
    Rd(i)=(und(i)-Rd(i))*deltaT/taod(i)+Rd(i-1);
while((sum(dd-Rd(i))>0.01)||(sum(dd-Rd(i))<-0.01))
         dd=Rd(i);
    und(i)=1/(1+exp(-(V(i)+10)/10));
    unf(i)=1/(1+exp((V(i)+48)/5));
    taod(i)=0.1*exp(-((V(i)+62)/13)^2)+0.05;
    taof(i)=0.5*exp(-((V(i)+55.6)/18)^2)+0.5;
    Rd(i)=(und(i)-Rd(i))*deltaT/taod(i)+Rd(i-1);
    Rf(i)=(unf(i)-Rf(i))*deltaT/taof(i)+Rf(i-1);
end;   
    JR(i)=gR*Rd(i).*Rf(i).*(V(i)-Eca(i));  %%%%JR is the Ca current via R type channel
end

%%% N-type %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function NTCC()   %%%%% T type channel, JT is the Ca current
global deltaT V i  gT  Eca Td ftf fts  JT dd 
if i==1
  Td(1)=1/(1+exp(-(V(1)+63.5)/1.5));
  ftf(1)=1/(1+exp((V(1)+76.2)/3));
  fts(1)=1/(1+exp((V(1)+76.2)/3));
  return;
else
  Td(i)=Td(i-1);
  ftf(i)=ftf(i-1);
  fts(i)=fts(i-1);
end;
und(i)=1/(1+exp(-(V(i)+63.5)/1.5));
ft(i)=1/(1+exp((V(i)+76.2)/3));
taodt(i)=65*exp(-((V(i)+68)/6)^2)+12;
taotf(i)=50*exp(-((V(i)+72)/10)^2)+10;
taots(i)=400*exp(-((V(i)+100)/10)^2)+400;
    Td(i)=(und(i)-Td(i))*deltaT/taodt(i)+Td(i-1);  
    dd=Td(i);
    Td(i)=(und(i)-Td(i))*deltaT/taodt(i)+Td(i-1);  
    ftf(i)=(ft(i)-ftf(i))*deltaT/taotf(i)+ftf(i-1);  
    fts(i)=(ft(i)-fts(i))*deltaT/taots(i)+fts(i-1);  
while((sum(dd-Td(i))>0.01)||(sum(dd-Td(i))<-0.01))
     dd=Td(i);
    und(i)=1/(1+exp(-(V(i)+63.5)/1.5));
   ft(i)=1/(1+exp((V(i)+76.2)/3));
   taodt(i)=65*exp(-((V(i)+68)/6)^2)+12;
    taotf(i)=50*exp(-((V(i)+72)/10)^2)+10;
   taots(i)=400*exp(-((V(i)+100)/10)^2)+400;
    Td(i)=(und(i)-Td(i))*deltaT/taodt(i)+Td(i-1);  
    ftf(i)=(ft(i)-ftf(i))*deltaT/taotf(i)+ftf(i-1);  
    fts(i)=(ft(i)-fts(i))*deltaT/taots(i)+fts(i-1);  
end;
    JT(i)=gT*Td(i).*(ftf(i)+0.04*fts(i)).*(V(i)-Eca(i));  %% JT is the Ca current via T type channel
end




