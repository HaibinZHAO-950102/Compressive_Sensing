# Compressive_Sensing  
Master Degree Thesis: Compressive Sensing + Filtering  

---------------------------------------------------
## 16.04.2020  
Today i wrote T_homo_modal_1.m, T_homo_modal_2.m, T_inhomo_modal_1.m, T_inhomo_modal_2.m to simulate the temperature distribution on a rod, using modal analysis.  

---------------------------------------------------
## 17.04.2020  
today i wrote a kalman filter to reconstructe the tempreature distribution on the rod.  
  
i met a problem, that i set the covariance-matirx of system-noise and measurement to 0. Thus the error-covariance-matrix will become zeros-matrix in the process of iteration, and then it is not invertible any more, which lead to the error of the programm.  
  
Tomorrow i will write another function to simulate the temperature distribution using fine difference method.  

------------------------------------------------------  
## 18.04.2020  
today i derived the formula in finite difference method, and then i wrote functions to simulate the temperature distribution using FDM in homogenous condition.  

------------------------------------------------------
## 19.04.2020
today i derived the formula of FDM in inhomogenous condition. and wrote functions to realize them. In addition, i wrote a function to compare 2 methods (FDM and MA).  
FDM need much longer time with the increasing number of nodes, but it works better than modal analysis around the points where the tempreature changes very fast.  
MA works faster but performance not good at the points where a dramatic change happens.  
  
people should pay attention to the delta-function in the numerical calculation. e.g. A delta function  
x = 3 * delta(t - 0.1) with t = [0, 0.1, 0.2, 0.3, 0.4]  
should be:  
x = [0, 3*(1/dt), 0, 0, 0]， in this example dt = 0.1, so that  
x = [0, 30, 0, 0, 0], instead of  
x = [0, 3, 0, 0, 0]  
  
Tomorrow I will study the two-dimensional distribution of temperature

----------------------------------------------------------
## 21.04.2020
today i derived the formula of Modal Analysis in 2D homogenous condition. It is more complicated and need more calculating time.  
then i wrote a function to realize it, and tried a example.
  
-------------------------------------------------------------
## 22.04.2020
today i derived the 2D inhomogenous situation using modal analysis. and wrote a function to realize it.  
in addition I tried 2 examples.  
tomorrow i will solve the same problem using FDM  

-------------------------------------------------------------
## 24.04.2020  
today i wrote some function to realize 2D tempreature distribution in homogenous condition. Using FDM  

-------------------------------------------------------------
## 25.04.2020  
today i realized 2D inhomogenous situation, using FDM  

-------------------------------------------------------------
## 26.04.2020  
today i begined to write some basic functions in Compressive Sensing, e.g. L1,L0,L2 Norm, something about sampling and so on.  

-------------------------------------------------------------
## 27.04.2020  
today i realized Compressive Sensing and it's Reconstruction for 1D Signal, tomorrow i will do some work for 2D Signal (e.g. for Image reconstruction)    

------------------------------------------------------------
## 29.04.2020
today i realized CS in 2D Situation (Image)  

------------------------------------------------------------
## 01.05.2020
today i improved my algorithm, using cosine transform as basic functions






