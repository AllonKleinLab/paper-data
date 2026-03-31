import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from tqdm import tqdm
import traceback
from pyDOE2 import lhs
import argparse


from typing import Optional, Callable, Dict, Any, Tuple
import types


from scipy.integrate import solve_ivp
from scipy.optimize import fsolve, minimize

from dataclasses import dataclass, field

import matplotlib.style as style
style.core.USER_LIBRARY_PATHS.append("/n/groups/klein/qiu/matplotlib_styles")

qw_light = '/n/groups/klein/qiu/matplotlib_styles/qw_light_paper.mplstyle'
plt.style.use(qw_light)

import warnings
warnings.filterwarnings("ignore")

### SAVING FILE NAME
# Set up argument parser
parser = argparse.ArgumentParser(description="Run simulation and save results with a custom file suffix.")
parser.add_argument("sim_num", type=str, help="Simulation identifier for the output filename")
args = parser.parse_args()
output_filename = f'simulation_results_114/simulation_tracker_df_5_{args.sim_num}.csv'
print(f"will save file: {output_filename}")

##### SUB ROUTINES

def get_gmin(A, N, tau):
    gmin = (1-A**(-1/N))*(N/tau)
    return gmin

# Function to run simulations for multiple parameter sets
# Function to run simulations for multiple parameter sets
def run_simulations(parameter_sets, t_span, t_eval):
    output_dynamics = {}
    for key, param_deviation in parameter_sets.items():
        simulation = Simulation()
        results = simulation.simulate(t_span, t_eval,param_deviation)
        output_dynamics[key] = results
        # print(f"E(0) in {key}: {results['E'][0]}")
    return output_dynamics

# This subroutine models loss of lung capacity - first gradual, then accelerating
def biphasic_alpha(instance, t):
    # A smooth function simulating loss of lung capacity over time
    # Can try different functional forms...
    return instance.alpha_min + (instance.alpha_max-instance.alpha_min)*np.exp(- (t/(4*instance.alpha_T))- 2*(t/(2*instance.alpha_T))**2)

# Function for mixture distribution for sampling parameters
def log_or_linear_random(min_value, max_value,p,loc, scale):
    lin = np.random.uniform(min_value,max_value)
    mid = np.abs(np.random.normal(loc=loc,scale=scale))
    log = np.exp(np.random.uniform(np.log(min_value+1e-2), np.log(max_value)))
    output = np.random.choice([lin,mid, log], p=p)
    return np.abs(output)

def pO2_performance_metric(output_dynamics,penalty=2,target=1):
    po2_t = np.array(output_dynamics['pO2'])
    t = np.array(output_dynamics['t'])
    performance_across_time = target - po2_t[:-1]
    dt = t[1:]-t[0:-1]
    #performance_across_time[performance_across_time < 0 ] = 0
    return np.sum(abs(performance_across_time)**(penalty) * t[:-1] * dt )**0.5
    
## Function for plotting output dynamics of simulations
def plot_variable_H2(output_dynamics,class_list):
    vars = ['pO2','x','y','surv','E']
    plt.figure(figsize=(10,2))
    for i,var in enumerate(vars):
        plt.subplot(1,5,i+1)
        for i,sim in enumerate(class_list):
            out = output_dynamics[sim]
            if i>=2:
                style = '--'
            else:
                style='-'
                
            if var == 'y':
                plt.plot(out['t'],out[var][-1],label=sim,linestyle=style,)
                #plt.yscale('log',base=10)
                #plt.ylim([0.0005,0.01])
            else:
                plt.plot(out['t'],out[var],label=sim,linestyle=style,)
                if var=='surv':
                    plt.ylim([0,1])
        plt.xlabel('Time (days)')
        plt.ylabel(var)
    plt.tight_layout()
    plt.legend(bbox_to_anchor=(1,1),fontsize=8)
    sns.despine()

## Function for plotting sensitivity simulations
def check_plots(tracker_df):
    normoxia = tracker_df[tracker_df.label=='param_testing_normoxia']
    hypoxia = tracker_df[tracker_df.label=='param_testing_hypoxia']
    plt.figure(figsize=(6,2))
    plt.subplot(1,3,1)
    sns.scatterplot(x=1/np.array(hypoxia.t_O2_rec),
                    y=1/(np.array(normoxia.progenitor_surv)+1e-2),
                    hue = np.array(normoxia.param_val),
                    palette='Reds',edgecolor='k',
              
            );
    
    plt.ylabel('Resting Cost\n(1/Freq. Prog. Surv.)')
    plt.xlabel('Speed:\n1/(Time to Target pO2 (Days))')
    plt.title(f'Tradeoff for adjusting {poi}', pad = 10,fontsize=8)
    plt.xlim([-0.01,0.15])
    plt.yscale('log')
    plt.legend([],[],frameon=False)
    sns.despine()
    
    plt.subplot(1,3,2)
    tracker_df = hypoxia.dropna()
    plt.plot(tracker_df.param_val,
             tracker_df.t_O2_rec,
             marker='o',markersize=5,
    
            );
    plt.xlabel(poi)
    plt.ylim([0,120])
    plt.ylabel('Time to Target pO2 (Days)')
    plt.title(f'Effect of {poi} on Recovery',pad=10,fontsize=8)
    sns.despine()
    
    plt.subplot(1,3,3)
    plt.plot(obtain_sensitivity_to_param(tracker_df,poi)[0],
             obtain_sensitivity_to_param(tracker_df,poi)[1],
             marker='o',
             markersize=5,
    
            );
    plt.xlabel(poi)
    plt.ylabel(f'dlog(t_pO2)/dlog({poi})')
    plt.title(f'Sensitivity',pad = 10,fontsize=8)
    sns.despine()
    plt.tight_layout()
    #plt.savefig(f'figures/{poi}_adjusted_values_of_interest.pdf')

    # Latin Hypercube Sampling to generate sample points within the parameter ranges
def latin_hypercube_sampling(param_range_dict, n_samples):
    param_names = list(param_range_dict.keys())
    n_params = len(param_names)
    
    # Generate LHS samples between 0 and 1
    lhs_samples = lhs(n_params, samples=n_samples)
    
    # Scale the samples to the actual parameter ranges
    param_samples = {}
    for i, param in enumerate(param_names):
        if param == 'K_r21':
            # Sample K_r21 in log space
            param_samples[param] = np.exp(np.log(param_range_dict[param][0]) + 
                                          (np.log(param_range_dict[param][1]) - np.log(param_range_dict[param][0])) * lhs_samples[:, i])
        else:
            # Sample other parameters in linear space
            param_samples[param] = param_range_dict[param][0] + \
                                   (param_range_dict[param][1] - param_range_dict[param][0]) * lhs_samples[:, i]
    
    return param_samples

# Function to check if normoxia dynamics are stable
def is_normoxia_stable(x_values, tolerance=1e-3):
    return np.abs(x_values[-1] - x_values[0]) < tolerance

# Function to calculate recovery speed for hypoxia condition
def calculate_recovery_speed(t_O2_rec):
    if pd.isna(t_O2_rec) or t_O2_rec == 0:
        return None
    return 1 / (t_O2_rec)


##### MODEL ####
@dataclass
class Simulation:
    N: float = 20 # Number of transit amplifying compartments
    
    epo_normalization: str = 'K1' # Set to either 'K1' or to 'E_0' to decide which will be set to 1.0
                                  # For scanning baseline designs, use 'E_0' (to keep baseline at E=1)
                                  # For physiological tuning, use 'K1' (because otherwise tuning will change K1)
    
    T: float = 40  # Lifetime of RBCs (days)
    tau: float = 7 # Mean time to propagate through progenitor pool

    K_r21: float = 2.0 #1.0 #3  # Ratio of K2/K1. Note that K2 is higher, because early compartment less sensitive

    gmin: float = 0.65 #0.8 #0.35 #0.7 #get_gmin(1000, 20, 7) # Baseline progenitor proliferaion at st-state
    del_g: float = 0.65 #1.1 #0.7 #0.5  #(gmax-gmin)
    #gmax: float = 2
    S_mpp = 1/8000
        
    # One of these two parameters is ignored and then fit in units of the other, (see 'epo_normalization')
    E_0: float = 1.0 ### Preset to 1 
    K1: float = 1.0    # Ks - Epo concentration to achieve 50% survival of final RBC progenitors
    
    n_s: int = 2
    n_a: int = 2 # hill coefficient of Epo response
   
    pSurv: float = 0.2 # Expected Epo-dependent survival rate of progenitors in steady-state
    # baselineAmp: float = 100 # Expected baseline amplification of progenitors

    
    # Trivial parameters:
    pO2_tgt: float = 1.0  # Target pO2
    fO2: float = 1.0      # External available O2
    xBar: float = 1.0 # Cell count at steady state (defined as 1.0)

    
    # Inflammation dynamics:
    alpha_max: float = 1.0 # Baseline absorbance (per unit x), pO2=fO2*x*alpha
    alpha_min: float = 1.0 # Minimum absorbance (lung infection) (per unit x) - when set to max, no disease
    alpha_T: float = 0.001  # Time (days) to reach minimum absorbance
    custom_alpha: Optional[Callable] = None  # Optional custom alpha function
                                            # alpha(self,t), to replace the default function
      
    # Feedback controller
    b: float = 5  # Sensitivity of Epo to pO2 deviation from target (dlog(Epo)/d(pO2))
    L_I: float = 0.05 # Strength of integral control
    integral_term: float = field(init=False, default=0.0)
    _last_t: float = 0.0

 
    # INITIALIZATION SUBROUTINES
    # ==========================
    def __init__(self):
            
        self.alpha_max = self.pO2_tgt / (self.fO2 * self.xBar)  # effective absorbance per RBC
        self.nu = self.N / self.tau  # 1/Lifetime in each transitional state (days)
        
    def __post_init__(self):
        self.alpha_max = self.pO2_tgt / (self.fO2 * self.xBar)  # effective absorbance per RBC
        self.nu = self.N / self.tau  # 1/Lifetime in each transitional state (days)
        
        # Replace alpha function with custom input if provided:
        if self.custom_alpha is not None:
            self.alpha = types.MethodType(self.custom_alpha, self)
        else:
            self.alpha = self.default_alpha

    # THE MODEL EQUATIONS
    # ==========================        
    def H(self, z, K, n):
        return np.abs(z) ** n / (np.abs(z) ** n + K ** n)

        
    def g(self, E_t, K2):
        #return self.gmin + (self.gmax - self.gmin) * self.H(E_t, K2, self.n_a)
        return self.gmin + self.del_g * self.H(E_t, K2, self.n_a)
    

    def S(self, E_t, K1, y10):
        return self.nu * self.H(E_t, K1, self.n_s) * y10


    def default_alpha(self, t):
        return self.alpha_min + (self.alpha_max - self.alpha_min) * np.exp(-t / self.alpha_T)

    def pO2(self, x, alpha_t):
        return self.fO2 * x * alpha_t

    def E(self, t, x):            
        pO2_diff = self.pO2(x, self.alpha(t)) - self.pO2_tgt
        if t == 0:
            self.integral_term = 0  # Reset integral term at t=0
        else:
            self.integral_term += pO2_diff * (t - getattr(self, '_last_t', 0))  # Accumulate the integral term
        self._last_t = t  # Update the last time point
        return self.E_0 * np.exp(-self.b * pO2_diff) - self.L_I * self.integral_term

    # The following are the dynamical equations themselves:
    def model(self, t, variables):
        x, *y = variables

        E_t = self.E(t, x)
        g_t = self.g(E_t, self.K2)

        dxdt = self.S(E_t, self.K1, y[-1]) - x / self.T

        dydt = [0] * self.N  # Initialize list of progenitor states
        dydt[0] = self.S_mpp - (self.nu - g_t) * y[0]  # First progenitor
        for i in range(1, self.N):
            dydt[i] = self.nu * y[i - 1] - (self.nu - g_t) * y[i]  # Subsequent progenitors

        return [dxdt] + dydt  # Convert into a list and return

    
    
    # SUBROUTINES TO FIT E0 OR K1 UNITS FROM STEADY-STATE REQUIREMENTS
    # ================================================================
    
    # IF FITTING E_0 WITH UNITS OF K1:
    def steady_state_equations_E0(self, E_0):
        g0 = self.g(E_0,self.K2)
        
        log_y = [np.log(self.S_mpp / (self.nu - g0))]
        for i in range(1, self.N):
            log_y.append( log_y[i - 1] + np.log(self.nu / (self.nu - g0)) )
        xBar = self.S(E_0, self.K1, np.exp(log_y[-1])) * self.T
        return (xBar - self.xBar)**2

    
    def calculate_E0(self):
        initial_guess = 1.0
        #E0_solution = fsolve(self.steady_state_equations, initial_guess,)[0]
        E0_solution = minimize(self.steady_state_equations_E0, 
                               x0=[initial_guess], bounds=[(1e-3, 20)], 
                               tol=1e-30, options={'maxiter': 1000}).x[0]
        return E0_solution

    
    # IF FITTING K1 WITH UNITS OF E_0:
    def steady_state_equations_K1(self, K1):
        K2 = self.K_r21 * K1
        g0 = self.g(self.E_0, K2)

        log_y = [np.log(self.S_mpp / (self.nu - g0))]
        for i in range(1, self.N):
            log_y.append( log_y[i - 1] + np.log(self.nu / (self.nu - g0)) )
        xBar = self.S(self.E_0, K1, np.exp(log_y[-1])) * self.T
        return (self.xBar - xBar)**2


    def calculate_K1(self):        
        initial_guess = 1.0
        K1_solution = minimize(self.steady_state_equations_K1, x0=[initial_guess], bounds=[(1e-3, 20)]).x[0]
        return K1_solution

    # THE SUBROUTINE THAT SETS THE MISSING PARAMTER (K1 OR E0) AND GETS STEADY-STATE:
    def initialize_steady_state(self):
        if self.epo_normalization == 'K1':
            self.K1 = self.calculate_K1()
            #print(f'{self.K1}')
            self.K2 = self.K_r21 * self.K1
        elif self.epo_normalization == 'E_0':
            self.K2 = self.K_r21 * self.K1
            self.E_0 = self.calculate_E0()
        else:
            raise ValueError("epo_normalization must be 'E_0' or 'K1' only.")
            
        
        g0 = self.g(self.E_0,self.K2)

        self.y_ss = [self.S_mpp / (self.nu - g0)]
        for i in range(1, self.N):
            self.y_ss.append(self.y_ss[i - 1] * self.nu / (self.nu - g0))
        
        # print(f"Steady state values: E_0 = {self.E_0}, xBar = {self.xBar}, y_ss = {self.y_ss}")
        return self.xBar, self.y_ss

    def solve_steady_state(self):
        g0 = self.g(self.E_0,self.K2)
  
        self.y_ss = [self.S_mpp / (self.nu - g0)]
        for i in range(1, self.N):
            self.y_ss.append(self.y_ss[i - 1] * self.nu / (self.nu - g0))
        
        return self.xBar, self.y_ss


    # POST-RUN DATA ANALYSIS SUBROUTINES
    # ==================================
    def find_time_for_midpoint(self, t, x):
        midpoint_value = 0.5 * (np.max(x) + np.min(x))

        for i in range(len(x) - 1, 0, -1):
            if (x[i] - midpoint_value) * (x[i - 1] - midpoint_value) <= 0:  # Crossing midpoint_value
                return t[i]
        return None  # Return None if no crossing is found

    def find_time_for_O2(self, t, x):
        target_value = 0.95
        
        for i in range(1, len(x)):
            if x[i - 1] < target_value <= x[i]:  # Crosses from below to above
                if x[i] == target_value:
                    return t[i]
                elif x[i - 1] == target_value:
                    return t[i - 1]
                else:
                    return t[i]
        return None  # Return None if no crossing is found  

    
    # THE SUBROUTINE THAT RUNS THE SIMULATION AND GENERATES POST-RUN ANALYSES
    # =======================================================================    
    def simulate(self, t_span, t_eval, perturbations: Optional[Dict[str, float]] = None):
        
        # Apply perturbations if provided
        if len(perturbations)>0:
            for key, value in perturbations.items():
                if hasattr(self, key):
                    setattr(self, key, value)
        
        self.__init__()

        # Initialize steady state
        self.xBar, self.y_ss = self.initialize_steady_state()
        #print('K1:',self.K1)
        
         # Ensure the model function uses the updated parameters
        initial_conditions = [self.xBar] + self.y_ss
        
        self.__post_init__()
        
        
        sol = solve_ivp(self.model, t_span, initial_conditions, t_eval=t_eval,  method='BDF')
        
        self._last_t = 0.0  # Reset integral term time tracking

        sol_x = sol.y[0]
        E_t = [self.E(t, x) for t, x in zip(sol.t, sol_x)]
        alpha_t = [self.alpha(t) for t in sol.t]
        pO2_t = [self.pO2(x, alpha) for x, alpha in zip(sol_x, alpha_t)]
        surv_t = [self.H(E, self.K1, self.n_s) for E in E_t]
        flux_t = [self.nu * y10 for y10 in sol.y[-1, :]]
        t_mid = self.find_time_for_midpoint(sol.t, pO2_t)
        t_O2 = self.find_time_for_O2(sol.t, pO2_t)

        self.results = {
            't': sol.t, 
            'x': sol_x, 
            'y': sol.y[1:], 
            'E': E_t, 
            'K1':self.K1,
            'alpha': alpha_t, 
            'pO2': pO2_t, 
            'surv': surv_t, 
            'flux': flux_t, 
            't_midO2': t_mid, 
            't_O2_rec': t_O2
        }
        
        return self.results



#### Simulation

# Initialize an empty list to collect simulation data and error logs
tracker_data = []
error_params = []


# Define the parameter range dictionary
param_range_dict ={
#'tau'  : [5,11], ## Mean time to propagate progenitor pool
'del_g':[0.25,1.5],
'K_r21' : [0.05,20], ## The ratio for which Epo surval over amplificaiton
'n_a': [1, 4], ## The hill function co-efficient that defines Epo reponse for survival
#'S_mpp': [(0.1/8000),(10/8000)], ## The flux of multipotent progenitors
'gmin': [0.3,1.0], ## The extent of Epo independent amplification
}

# Target number of valid simulations
target_valid_simulations =1000
valid_simulations = 0
total_simulations_run = 0  # Track total simulations attempted

# Create a progress bar for valid simulations
with tqdm(total=target_valid_simulations, desc="Valid Simulations") as pbar:
    while valid_simulations < target_valid_simulations:
        # Generate parameter samples in batches
        param_samples = latin_hypercube_sampling(param_range_dict, target_valid_simulations - valid_simulations)

        # Run simulations for the sampled parameter sets
        for i in range(len(param_samples['K_r21'])):
            total_simulations_run += 1
            param_values = {key: param_samples[key][i] for key in param_range_dict}

            parameter_sets = {
                'normoxia_il17': {
                    'gmin': param_values['gmin'],
                    'del_g': param_values['del_g'],
                    'K_r21': param_values['K_r21'],
                    'n_a': param_values['n_a'],
                    #'S_mpp': param_values['S_mpp'],
                    #'tau': param_values['tau']
                },
                'hypoxia_il17': {
                    'gmin': param_values['gmin'],
                    'del_g': param_values['del_g'],
                    'K_r21': param_values['K_r21'],
                    'n_a': param_values['n_a'],
                    #'S_mpp': param_values['S_mpp'],
                    #'tau': param_values['tau'],
                    'alpha_min': 0.7, 
                    'alpha_T': 5, 
                    'custom_alpha': biphasic_alpha
                }
            }
            
            t_span = (0, 100)
            t_eval = np.linspace(t_span[0], t_span[1], 2000)
            
            try:
                
                # Run the simulations for both normoxia and hypoxia conditions
                output_dynamics = run_simulations(parameter_sets, t_span, t_eval)
                
                # Check stability for normoxia condition
                normoxia_stable = is_normoxia_stable(output_dynamics['normoxia_il17']['x'])
                valid_hypoxia = True
                temp_tracker_data = []
                
                for condition in ['normoxia_il17', 'hypoxia_il17']:
                    flux = output_dynamics[condition]['flux'][0]
                    surv = output_dynamics[condition]['surv'][0]
                    
                    if surv <= 0 or surv > 1:
                        print('Survival is in an unreal range')
                        valid_hypoxia = False
                        break
                        
            

                    peak_flux = np.max(output_dynamics[condition]['flux'])
                    min_pO2 = np.min(output_dynamics[condition]['pO2'])
                    K1 = np.min(output_dynamics[condition]['pO2'])
                    mid_pO2_t = output_dynamics[condition]['t_midO2']
                    t_O2_rec = output_dynamics[condition]['t_O2_rec']
                    recovery_speed = calculate_recovery_speed(t_O2_rec)
                    resting_cost = 1/output_dynamics[condition]['surv'][0]
                    penalty_pO2_t = pO2_performance_metric(output_dynamics[condition], target=1)

                    if condition == 'normoxia_il17' and resting_cost>40 :
                        #print(f'Resting Cost is too high: {resting_cost}')
                        valid_hypoxia = False
                        break
                    elif resting_cost is None:
                        #print('there is no valid resting cost')
                        valid_hypoxia = False
                        break
                    #print('I finished here')
                    temp_tracker_data.append({
                        'simulation_number': total_simulations_run,
                        'del_g': param_values['del_g'],
                        'K_r21': param_values['K_r21'],
                        'gmin': param_values['gmin'],
                        #'tau': param_values['tau'],
                        'n_a': param_values['n_a'],
                        #'S_mpp':param_values['S_mpp'],
                        'K1': output_dynamics[condition]['K1'],
                        'alpha_min': parameter_sets[condition].get('alpha_min', np.nan),
                        'alpha_T': parameter_sets[condition].get('alpha_T', np.nan),
                        'progenitor_flux': flux,
                        'progenitor_surv': surv,
                        'peak_flux': peak_flux,
                        'min_pO2': min_pO2,
                        'mid_pO2_t': mid_pO2_t,
                        't_O2_rec': t_O2_rec,
                        'recovery_speed': recovery_speed,
                        'resting_cost': resting_cost,
                        'condition': condition,
                        'normoxia_stable': normoxia_stable if condition == 'normoxia_il17' else np.nan,
                        'penalty_pO2_t': penalty_pO2_t
                    })
                #print('I finished appending')
                # Add data if both conditions are valid
                if valid_hypoxia and normoxia_stable:
                    
                    tracker_data.extend(temp_tracker_data)
                    valid_simulations += 1
                    pbar.update(1)
            
            except Exception as e:
                continue

            # Stop if the target number of valid simulations is reached
            if valid_simulations >= target_valid_simulations:
                break
                      
###

# Convert the collected data into a DataFrame
tracker_df = pd.DataFrame(tracker_data)

# Save the tracker_df to file
tracker_df.to_csv(output_filename, index=False)

print(f"Total simulations run: {total_simulations_run}")
print(f"Valid simulations retained: {valid_simulations}")