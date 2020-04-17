# Compressive_Sensing
Master Degree Thesis: Compressive Sensing + Filtering  

---------------------------------------------------
## 16.04.2020
Today i wrote T_homo_modal_1.m, T_homo_modal_2.m, T_inhomo_modal_1.m, T_inhomo_modal_2.m to simulate the temperature distribution on a rod, using momodal analysis.  

---------------------------------------------------
## 17.04.2020
today i wrote a kalman filter to reconstructe the tempreature distribution on the rod. i met a problem, that a set the covariance-matirx of system-noise and measurement to 0. Thus the error-covariance-matrix will become zeros-matrix in the process of calculation, and then it is not invertible any more, which lead to the error of the programm.  
Tomorrow i will write another function to simulate the temperature distribution using fine difference method.  

------------------------------------------------------
