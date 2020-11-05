# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 18:39:47 2020

@author: Ariadna
"""

#%% importamos librerias
import numpy as np
import scipy.stats as si
import sympy as sy
from sympy.stats import Normal, cdf
#%%
So=100
K=100
r=0.02
T=1
sigma=0.3
#función PUT
def E_put(So,K,r,sigma,T):    
    d1=(np.log(So/K)+(r+0.5*sigma**2)*T)/(sigma*np.sqrt(T))
    d2=(np.log(So/K)+(r-0.5*sigma**2)*T)/(sigma*np.sqrt(T))

    Put=K*np.exp(-r*T)*si.norm.cdf(-d2)-So*si.norm.cdf(-d1)
    
    return Put

print( E_put(So,K,r,sigma,T))
#función CALL
def E_call(So,K,r,sigma,T):    
    d1=(np.log(So/K)+(r+0.5*sigma**2)*T)/(sigma*np.sqrt(T))
    d2=(np.log(So/K)+(r-0.5*sigma**2)*T)/(sigma*np.sqrt(T))

    Call=So*si.norm.cdf(d1)-K*np.exp(-r*T)*si.norm.cdf(d2)
    
    return Call

print(E_call(So,K,r,sigma,T))
#%% delta
def Delta_PUT(So,K,r,sigma,T):
    d1=(np.log(So/K)+(r+0.5*sigma**2)*T)/(sigma*np.sqrt(T))
    delta=-si.norm.cdf(-d1)
    return delta

print(Delta_PUT(So,K,r,sigma,T))
#%% gamma
def Gamma_PUT(So,K,r,sigma,T):
    d1_negativo=-(np.log(So/K)+(r+0.5*sigma**2)*T)/(sigma*np.sqrt(T))
    N_d1_negativo=(1/np.sqrt(2*np.pi))*np.exp(-.5*d1_negativo**2)
    gamma=N_d1_negativo/(So*sigma*np.sqrt(T))
    return gamma

print(Gamma_PUT(So,K,r,sigma,T))
#%% rho
def Rho_PUT(So,K,r,sigma,T):
    d2=(np.log(So/K)+(r-0.5*sigma**2)*T)/(sigma*np.sqrt(T))
    rho=-T*K*np.exp(-r*T)*si.norm.cdf(-d2)
    return rho

print(Rho_PUT(So,K,r,sigma,T))
#%% theta
def Theta_PUT(So,K,r,sigma,T):
    d1_negativo=-(np.log(So/K)+(r+0.5*sigma**2)*T)/(sigma*np.sqrt(T))
    d2=(np.log(So/K)+(r-0.5*sigma**2)*T)/(sigma*np.sqrt(T))
    N_d1_negativo=(1/np.sqrt(2*np.pi))*np.exp(-.5*d1_negativo**2)
    theta=r*K*np.exp(-r*T)*si.norm.cdf(-d2)-.5*(sigma/np.sqrt(T))*So*N_d1_negativo
    return theta

print(Theta_PUT(So,K,r,sigma,T)) 
#%% vega
def Vega_PUT(So,K,r,sigma,T):
    d1=(np.log(So/K)+(r+0.5*sigma**2)*T)/(sigma*np.sqrt(T))
    N_d1=(1/np.sqrt(2*np.pi))*np.exp(-.5*d1**2)
    vega=So*N_d1*np.sqrt(T)
    return vega
print(Vega_PUT(So,K,r,sigma,T))