##intermediate variables  

# Economy
kapa    = max(kmin,min(kmax,k_0+k_pi*pi))                                       # Investment function
I       = multiply*kapa*Y
div     = max(div_min, min(div_max, div_0+div_pi*pi_k))                         # Dividend function
phi     = phi_0+phi_1*lambda                                                    # Philips Curve
beta    = q * (1 - N/N_max)                                                     # Population growth rate
L       = Y_0/a                                                                 # Labor
lambda  = L/N                                                                   # Employment
omega   = W*L/(p*Y)                                                             # Wage share 
c       = W*L/(p*Y)                                                             # Production cost

delta_k = delta + dam_k                                                         # Damage and depreciation on capital

# Output
Y_0     = K/v                                                                   # Output net of damage & abatement
Y       = Y_0*(1-dam_y)*(1-A)-MC

# Debt & Profit
d       = D/(p*Y)                                                               # debt to output ratio
i       = eta*(mu*(c+omitted)-1)                                                # Inflation 
PI      = p*Y-W*L-r*D-p*delta_k*K-p*T_gc                                        # Profit
pi      = max(0,PI/(p*Y))                                                              # Profit share
pi_k    = PI/(p*K)                                                              
Div     = div*p*K                                                               # Dividends
PI_r    = PI-Div                                                                # Profit net of dividends
 
# Abatement 

A       = (sigma*p_bs/(1000*theta))*n^theta                                     # Abatement
n       = min((p_car/p_bs)^(1/(theta-1)),1)                                     # Emission reduction rate
T_gc    = tax*conv10to15*FO_c                                                   # Carbon tax on fossil fuel

# Climate
dam     = dam_on*(1-1/(1+pi_1*Temp+pi_2*Temp^2+pi_3_weitz*Temp^exp_NQ))         # Damage (dam_on: 1 = damages active, 0 = no-damage scenario)
dam_k   = f_K*dam                                                               # Capital damage (f_K: 0 = low damage / output only, 1/3 = high damage / output+capital)
dam_y   = 1-(1-dam)/(1-dam_k)                                                   # Production damage

#E_ind   = Y_0*sigma*(1-n)                                                      # Industrial emission
#E_tro   = lambda_1*(1-lambda_2)^t                                              # defined by a time derivative Carbon emission from tropical forest
EE      = FO_c + HB_bor*theta_bor+HB_tem*theta_tem+HB_tro*theta_tro             # Emission from energy conversion 
ES = CH*delta_CH+(1-v_c)*(HS_bor*theta_bor+HS_tem*theta_tem+HS_tro*theta_tro)   # Release of carbon from forest products
EF_bor  = -theta_bor*(phi_bor*F_bor*(1-F_bor/Fmax_bor))                         # CO2 Forest sequestration (boreal)
EF_tem  = -theta_tem*(phi_tem*F_tem*(1-F_tem/Fmax_tem))
EF_tro  = -theta_tro*(phi_tro*F_tro*(1-F_tro/Fmax_tro)-D_tro)
EF      = EF_bor+EF_tem+EF_tro                                                  # Total forest sequestration
E       = EE + ES + EF                                                          # Total carbon emission
E_ant   = E-EF	
	
F_ind   = (F2xCO2/log(2))*log(CO2_at/CAT_pind)                                  # Industrial forcing
F_exo   = 0.5 + t/(time_sim*2)                                                  # Exogenous radiative forcing ; Linear ; [0.5;1]
Forcing = F_ind + F_exo                                                         # Radiative forcing

C       = 5/(1/C_init+beta_C*(S-S_ref))                                         # Adaptative function of climate
tax     = min(p_car,p_bs)   

# Energy & emission

Energy   = ksi*Y_0*sigma*(1-n)                                                  # Energy to produce Y
FO_c     = min(((p_FO+p_car)/p_BI *(1-beta_car)/beta_car)^(beta_car-1)*Energy/zeta,FO_r) # Fossil fuel demand from Cobb-Douglas
                                            

F_tro_per = 100*F_tro/388

BI_m3    = (Energy/zeta)^(1/(1-beta_car))*FO_c^(-beta_car/(1-beta_car))+max(FO_c-FO_r,0)     # Bioenergy demand

HB_bor   = BI_m3 * PHI                                                          # Boreal biomass harvested for bioenergy 
HB_tem   = BI_m3 * (1-k-PHI)                                                    # Temperate biomass harvested for bioenergy
HB_tro   = BI_m3 * k  

# Forest
# Deforestation control rate RD (= DC in the paper), driven from the policy dashboard in run_model.R:
#   RD_slope = 0     -> Baseline (no forest policy)
#   RD_slope = 1/84  -> RD reaches 1 (halts net deforestation) in 2100
#   RD_slope = 1/42  -> RD reaches 2 in 2100  (Forest Policy A)
#   RD_slope = 1/22  -> RD reaches ~3.8 in 2100 (Forest Policy B, stronger reforestation)
RD = RD_slope*t                                                                # Deforestation control rate
D_tro = (1-RD)*E_tro/theta_tro                                                  # Total biomass deforestation
RE = E_tro*RD     	                                                            # Reduction of direct carbon emission from deforestation


H_bor = HS_bor + HB_bor                                                         # Total boreal harvesting
H_tem = HS_tem + HB_tem
H_tro = HS_tro + HB_tro

HS = HS_bor + HS_tem + HS_tro                                                   # Total biomass harvested for forest products

delta_CH = log(2)/HL                                                             # Fraction of carbon lost from the carboon stock in wood product per year
Inflow = v_c*(HS_bor*theta_bor+HS_tem*theta_tem+HS_tro*theta_tro)               # Inflow of carbon in wood products

# Marginal vs total forest-policy cost.
# f1,f2 are the Kindermann marginal-abatement-cost fit: MC_marg(RE)=f1*exp(f2*RE) ($/tC).
# The total cost entering GDP is its integral, C_F(RE)=(f1/f2)*(exp(f2*RE)-1), so that it
# vanishes at RE=0 and its slope reproduces the fitted Kindermann curve exactly.
#   MC_on   = 1 -> policy cost active ; 0 -> C_F = 0 (baseline)
#   MC_mode = 1 -> integral/total form C_F=(f1/f2)*(exp(f2*RE)-1) (zero at RE=0)
#   MC_mode = 0 -> level form (f1/f2)*exp(f2*RE) (diagnostic only)
MC =  MC_on*(f1/f2)*(exp(f2*RE)-MC_mode)
#MC = f1*RE^f2+((f3+f4*t)^f5*RE-1)                                               # Marginal cost function (do not deforest)

##time derivatives
a        = alpha*a                                                              # Productivity 
W        = phi*W                                                                # Wages
D        = p*I-PI_r-(p*delta_k*K)                                              # Debt
N        = N*beta                                                               # Workforce
K        = I-delta_k*K                                                       # Capital

sigma    = g_sigma*sigma                                                        # Growth of emission intensity
g_sigma  = delta_g_sigma*g_sigma

p        = i*p                                                                  # Price 
p_bs     = g_pbs*p_bs                                                           # Backstop technology price
p_car    = g_p_car*p_car                                                        # Carbon price

# CO2 Emissions
E_tro   = delta_E_tro*E_tro                                                     # Emission from deforestation 

CO2_at   = E/3.666-CO2_at*phi_12+phi_12*CO2_up*CAT_pind/CUP_pind
CO2_up   = phi_12*CO2_at-(phi_12*CAT_pind/CUP_pind+phi_23)*CO2_up+phi_23*CO2_lo*CUP_pind/CLO_pind
CO2_lo   = phi_23*CO2_up-phi_23*CO2_lo*CUP_pind/CLO_pind

# Temperature
Temp      = Forcing/C - F2xCO2/(C*S)*Temp - (gamma_ast/C)*(Temp-Temp_lo) 
#Temp     = (Forcing - Temp*F2xCO2/S-gamma_ast*(Temp-Temp_lo))/C                # Atmosphere and upper ocean mean temperature
Temp_lo  =  gamma_ast*(Temp-Temp_lo)/C_0                                        # Deep ocean mean temperature

# Forest

#Fmax_bor = -Fmax_bor*D_bor/F_bor                                                # Forest biomass carrying capacity
#Fmax_tem = -Fmax_tem*D_tem/F_tem
Fmax_tro = -Fmax_tro*D_tro/F_tro

F_bor = phi_bor*F_bor*(1-F_bor/Fmax_bor)-H_bor                                  # Forest biomass dynamics
F_tem = phi_tem*F_tem*(1-F_tem/Fmax_tem)-H_tem
F_tro = phi_tro*F_tro*(1-F_tro/Fmax_tro)-H_tro-D_tro

HS_bor = x_bor*HS*(((I-delta_k*K-MC)*a-K*alpha*a)/(K*a)+1)-HS_bor                           # Harvested boreal biomass for wood products dynamics
HS_tem = x_tem*HS*(((I-delta_k*K-MC)*a-K*alpha*a)/(K*a)+1)-HS_tem
HS_tro = x_tro*HS*(((I-delta_k*K-MC)*a-K*alpha*a)/(K*a)+1)-HS_tro

CH = -delta_CH*CH+Inflow
# HS = HS*((I-delta_k*K-TC)/(L*v*a)+1)

FO_r = -FO_c

#TC = (1/10)*(f1/(f2+1))*RE^(f2+1)+1/(f5+log(f3+f4*(t+1)))*((f3+f4*(t+1))^(f5*RE)-1)-RE  # Total cost of avoiding deforestation


##initial values

a = 18.377
W = 10.608                                                                      # 0.578*59.7/3.26025 with 3.26 = L = lambda*N
D  = 91.569                                                                     # debt Y*1.53
N = 4.8255                                                                      # Workforce (bill)
K = 161.57
# Capital : v*Y/((1-dam_Y_0)*(1-A_O))
sigma = 0.61867                                                                 # E_ind/((1-n)*Y) current emission intensity of the economy
g_sigma = - 0.0105                                                              # Growth rate of the emission intensity of the economy
p = 1
p_bs = 547.22
p_car = 2.0024
E_tro = 3.3                                                                     # Exogenous land use change CO2-e emissions http://www.fao.org/faostat/en/#data/GF
CO2_at = 851                                                                    # CO2e concentration in the atmosphere layer
CO2_up = 460                                                                    # CO2e concentration in the biosphere and upper ocean layer
CO2_lo = 1740                                                                   # CO2e concentration in the deeper ocean layer
Temp = 0.85                                                                     # Atmosphere, upper ocean layer, biosphere temperature
Temp_lo = 0.0068                                                                # Deeper ocean layer temperature


Fmax_tro = 776

F_bor = 173
F_tem = 84
F_tro = 388

HS_bor = 0.7649
HS_tem = 0.6599
HS_tro = 0.4252

CH = 4.8                                                                        # Stock of carbon in wood product 2015
FO_r = 7000 
#MC = 0

##parameters
# Fundamental Constants
alpha = 0.02                                                                    # Constant growth rate of productivity
a_0 = 18                                                                        # Productivity
q = 0.0305                                                                      # Constant growth rate of workforce 
q_g = 0.02744369                                                                # Constant growth rate of population
delta = 0.04                                                                    # Depreciation rate of capital
v = 2.7                                                                         # Constant Capital to Output ratio
r = 0.03                                                                        # Short term interest rate of the economy
g_pbs = -.00505076337                                                           # Exogenous growth rate of the back-stop technology price
NG_max = 12                                                                     # Maximum population (bill)
N_max = 7.055925493
a_pc = 0.                                                                       # a and b : parameters of carbon price trajectory 
b_pc = 0.02                                                                     # carbon-price trajectory kept below 2.5 °C in 2100
g_p_car = scenario_g_p_car                                                      # Carbon-price growth rate (set by the policy dashboard in run_model.R)
time_sim = 84                                                                   # 2016-2100 
conv10to15 = 1.160723971/1000

# ---- Scenario switches (all set from the POLICY DASHBOARD in run_model.R) ----
RD_slope = scenario_RD_slope                                                    # Slope of the deforestation-control rate RD = RD_slope*t (0 = Baseline)
MC_on    = scenario_MC_on                                                       # 1 = forest-policy cost active, 0 = no forest policy (C_F = 0)
MC_mode  = scenario_MC_mode                                                     # 0 = level C_F ; 1 = incremental C_F (zero at t0)
dam_on   = scenario_dam_on                                                      # 1 = climate damages active, 0 = no-damage scenario

# Investment Function
k_0 = .031774                                                                   # Constant of the investment function
k_pi = .575318414                                                               # Slope of the investment function
kmin = 0                                                                        # Lower bound of the investment function
kmax = 0.3                                                                      # Upper bound of the investment function
multiply = 1/.8                                                                 # Calibration scale on investment: pins the model's I/Y to the empirical GFCF/GDP ratio (0.2174) at the 2015 initial condition

# Philips Curve 
phi_0 = - 0.291535421                                                           # Constant of Phillips Curve function
phi_1 = 0.468777035                                                             # Slope of Phillips Curve function

#Inflation 
eta = 0.3                                        
mu = 1.2                                                                        # Relaxation parameter of the inflation
omitted = 0.3                                                                   # Flat proxy for the non-wage components of Bovari (2018) unit cost (Eq. 11): here c = wage share only, so 'omitted' stands in for capital + carbon unit costs
#conv10to15 = 1.160723971/1000                                                  # 1.160723971 conversion from 2010->2015 and /1000 scale factor with trillon dollars


# Dividend Function 
div_0 = .0512863
div_pi = .4729
div_min = 0
div_max = 1

# Climate 

# Climate Damages
pi_1 = 0
pi_2 = .00236
pi_3_weitz = .00000507                                                          # 0.000005070
pi_3_stern = .0000819
exp_NQ = 6.754
f_K = scenario_f_K                                                               # Fraction of environmental damage allocated to capital: 0 = low damage, 1/3 = high damage (set by dashboard)
theta = 2.6                                                                     # Parameter of the abatement cost function

# Temperature 
C_init = 1/0.098                                                                # Heat capacity of the atmosphere, biosphere and upper ocean
C_0 = 3.52                                                                      # Heat capacity of the deeper ocean
T_pind = 13.74                                                                  # Preindustrial temperature
S = 3.1                                                                         # Equilibrium climate sensitivity
S_ref = 2.9
beta_C = .01243
gamma_ast = 0.0176                                                              # Heat exchange coefficient between temperature layers

# CO2 parameters 
CAT_pind = 588                                                                  # CO2-e preindustrial concentration in the atmosphere layer
CUP_pind = 360                                                                  # CO2-e preindustrial concentration in the biosphere and upper ocean layer
CLO_pind = 1720                                                                 # CO2-e preindustrial concentration in the deeper ocean layer

phi_12 = .0239069                                                               # Transfer coefficient for carbon from the atmosphere to the upper ocean/biosphere
phi_23 = .0013409                                                               # Transfer coefficient for carbon from the upper ocean/biosphere to the lower ocean

delta_g_sigma = -.001                                                           # Variation rate of the growth of emission intensity 
delta_E_tro = -.0220096                                                         # Growth rate of deforestation change CO2 emissions

ksi = 14.39687278                                                               # Carbon emission energy parameter 
	
HL = 40                                                                         # Half life years for carbon in wood products
delta_CH = 0.0173287                                                            # Fraction of carbon lost from the carboon stock in wood product per year
car_wood = 0.18                                                                 # Share of carbon in roundwood harvest to wood products

# RadiativeForcing 

F2xCO2 = 3.6813                                                                  # Change in the radiative forcing resulting from a doubling of CO2-e concentration w.r.t. to the pre-industrial period
F_exo_start = 0.5                                                               # Initial value of the exogenous radiative forcing
F_exo_end = 1                                                                   # Value of the exogenous radiative forcing in 2100
F_exo_end_year = 2100                     
F_exo_end_start= 2016


# Marginal cost

f1 = 1.52638276079                                                                      # Cost parameter of reducing emissions from avoided deforestation
f2 = 1.18258619118984

elast = 2.0                                                                     # Elasticity of marginal utility of consumption
rho = 0.015                                                                     # Pure rate of social time preference

# Forest 

Fmax_bor = 346
Fmax_tem = 168

phi_bor = 0.01379341861                                                         # Intrinsic growth rate of boreal forest biomass
phi_tem = 0.03362579279                                                         # From forest dice, adapted to a continuous
phi_tro = 0.04051984059                                                         # framework 

theta_bor = 0.409                                                               # Carbon intensity in boreal forest biomass
theta_tem = 0.524
theta_tro = 0.611

x_bor = 0.41                                                                    # boreal biomass share of total roundwood harvest
x_tem = 0.36
x_tro = 0.23


HL = 40                                                                         # Half-life years for carbon in wood products
#HS = 1.85                                                                       # Harvested biomass for products 2015

v_c = 0.18                                                                      # Share of carbon in roundwood harvest to wood products

# Energy

beta_car = 0.902                                                                # Fossil fuel carbon elasticity 
zeta =20.72397771                                                          # Energy parameter


p_FO = 34.024
p_BI = 64.20567885

F0_max = 7000                                                                   # Maximum consumption fossil fuels (billions of tc)
EtoC  = 0.018926
Etom3 = 0.043939
w_bio = 1.902                                                                   # Bioenergy parameter
PHI = 0.068                                                                     # Bioenergy harvest elasticity boreal forest biomass
k = 0.7917                                                                      # Bioenergy harvest elasticity tropical forest biomass 




#_______________________________________________________________________________________________________________________________#
##time
begin = 0
end = 84
by = 0.05