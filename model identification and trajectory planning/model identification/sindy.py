#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Data:22/02/2021
Version:1.2
Author: Minyu Zhang, Jiangqiu Chen, Zhifei Yang
Project: Obtacle avoidancea and energy optimization for 3DOF(RRR) manipulator
Function of python code:
    Measurement noise and denoise simulation.
    Sparse model regression based on measurement data.
    Auto alpha selection for L1 regulation.
    Correlation filter reduce library size.
Code Environments:
    numpy==1.19.5
    scikit-learn==0.24.1
    scipy==1.5.4
    matplotlib==3.3.2
"""

import numpy as np
from sklearn import linear_model
from sklearn.metrics import mean_absolute_percentage_error as mape
from sklearn.metrics import mean_squared_error as mse
from noise import *
from _decimal import Decimal
import scipy.io as sio
import matplotlib.pyplot as plt
np.set_printoptions(linewidth=1000) # maximum print width in one signle line

cindex = { # color index
'aliceblue':            '#F0F8FF',
'antiquewhite':         '#FAEBD7',
'aqua':                 '#00FFFF',
'aquamarine':           '#7FFFD4',
'azure':                '#F0FFFF',
'beige':                '#F5F5DC',
'bisque':               '#FFE4C4',
'black':                '#000000',
'blanchedalmond':       '#FFEBCD',
'blue':                 '#0000FF',
'blueviolet':           '#8A2BE2',
'brown':                '#A52A2A',
'burlywood':            '#DEB887',
'cadetblue':            '#5F9EA0',
'chartreuse':           '#7FFF00',
'chocolate':            '#D2691E',
'coral':                '#FF7F50',
'cornflowerblue':       '#6495ED',
'cornsilk':             '#FFF8DC',
'crimson':              '#DC143C',
'cyan':                 '#00FFFF',
'darkblue':             '#00008B'}
# color dictionary for figure plotting
color=[]
for i in cindex.keys():
    color.append(i)
# ==========================================================================================================================================
"""
Training data generation
"lib" stores numerical training data
"libname" stores coresponding features
"""
X=sio.loadmat('data.mat')['x'][:,3:12]
p=sio.loadmat('data.mat')['x'][:,0:3]
##t = np.linspace(0, 2, 10001)
lib=[]
libname=[]

a=np.array([X[:,6],X[:,7],X[:,8],X[:,3]*X[:,4],X[:,4]*X[:,5],X[:,3]*X[:,5]])
a=a.T
b=np.array([X[:,3]*X[:,3],X[:,4]*X[:,4],X[:,5]*X[:,5]])
b=b.T
c=np.array([np.cos(X[:,1]),np.cos(X[:,2]),np.cos(2*X[:,1]),np.cos(2*X[:,2]),np.cos(X[:,1]+X[:,2]),np.cos(2*X[:,1]+X[:,2]),np.cos(2*X[:,1]+2*X[:,2])])
c=c.T
d=np.array([np.sin(X[:,2]),np.sin(2*X[:,1]+X[:,2]),np.sin(2*X[:,1]+2*X[:,2])])
d=d.T
e=np.array([np.sin(2*X[:,2])])
e=e.T

aname=['y7','y8','y9','y4*y5','y4*y6','y5*y6']
bname=['y4^2','y5^2','y6^2']
cname=['cos(y2)','cos(y3)','cos(2*y2)','cos(2*y3)','cos(y2+y3)','cos(2*y2+y3)','cos(2*y2+2*y3)']
dname=['sin(y3)','sin(2*y2+y3)','sin(2*y2+2*y3)']
ename=['sin(2*y2)']      
# ===========================================================================================================================================

for i in range(a.shape[1]):
    for j in range(c.shape[1]):
        lib.append(a[:,i]*c[:,j])
        libname.append(aname[i]+'*'+cname[j])

for i in range(b.shape[1]):
    for j in range(d.shape[1]):
        lib.append(b[:,i]*d[:,j])
        libname.append(bname[i]+'*'+dname[j])

for i in range(a.shape[1]):
    for j in range(d.shape[1]):
        lib.append(a[:,i]*d[:,j])
        libname.append(aname[i]+'*'+dname[j])

lib.append(c[:,0])
libname.append(cname[0])
lib.append(c[:,4])
libname.append(cname[4])
lib.append(b[:,0]*e[:,0])
libname.append(bname[0]+'*'+ename[0])

for i in range(3):
    lib.append(a[:,i])
    libname.append(aname[i])
### =====include the constant feature or no =====  
##lib.append(np.ones([len(X)]))
##libname.append('1')
### =============================================    
lib1=lib.copy()
lib2=lib.copy()
lib3=lib.copy()
lib1=np.array(lib1).T
lib2=np.array(lib2).T
lib3=np.array(lib3).T
# =============================================================================================================================================


"""
Measuremnet Nosie simulation & denoise algorithm
"pn" is the simulated noisy data
"NoiseEs" is estimation of measuremnet noise
"pes" is estimation of training data after denoise
"""
NoiseLevel=5
np.random.seed(0)

# Generate the noise
NoiseMag=[np.std(p[:,i])*NoiseLevel*0.01 for i in range(p.shape[1])]
Noise=np.hstack([NoiseMag[i]*np.random.randn(p.shape[0],1) for i in range(p.shape[1])])
# Add the noise and get the noisy data
pn=p+Noise

NoiseEs,pes=approximate_noise(np.transpose(pn), 1)
pes=np.transpose(pes)
mse_error=mse(p,pn)
mse_denoise=mse(p,pes)
p_mean=np.mean(p)
error_mean=np.mean(pn)
error_est=np.mean(pes)
p_std=np.std(p)
error_std=np.std(pn)
error_est_std=np.std(pes)
print('mse for error is %.10f'%(mse_error))
print('mse for error_est is %.10f'%(mse_denoise))
print(p_mean,p_std)
print(error_mean,error_std)
print(error_est,error_est_std)

# figure for visulization of denoise algorithm
plt.figure(1)
plt.plot(pn[:,0],linewidth=0.5,color='r',label='noise')
plt.plot(pes[:,0],linewidth=0.8,color='k',label='denoisy') ##plot denoisy data
plt.title('noise level=5%')
plt.xlabel("sample indice")
plt.ylabel("τ1(N/m)")
plt.legend(['measurement','denoisy'])
plt.grid()
# legend matrix for plotting
t1_legend=['denoisy']
t2_legend=['denoisy']
t3_legend=['denoisy']
# =============================================================================================================================================
##
##"""
##Main iterations include
##"iter_alpha" times traversals to looking for the best alpha for L1 regulation
##"turns" times iterations sparsing libraries
##"""
##iter_alpha=5 
##turns=30
##mapelist=np.zeros([3,iter_alpha])
##alphalist=np.zeros([3,iter_alpha])
##for al in range(iter_alpha+1):   # include one more round for best result output
####    alpha1=0.000005
####    alpha2=0.00005
####    alpha3=0.00004
##    if al in range(iter_alpha):
##        alpha1=float(Decimal(str(0.00000025))*Decimal(str((al+1))))
##        alpha2=float(Decimal(str(0.00001))*Decimal(str((al+5))))
##        alpha3=float(Decimal(str(0.00001))*Decimal(str((al+4))))
##    else:
##        alpha1=bestalpha[0]
##        alpha2=bestalpha[1]
##        alpha3=bestalpha[2]
##        print('\r\n')
##        print('best alpha is(α1,α2,α3)=(%8f,%8f,%8f)'%(alpha1,alpha2,alpha3))
##        print('\r\n')
##    reg1=linear_model.Lasso(alpha=alpha1)
##    reg2=linear_model.Lasso(alpha=alpha2)
##    reg3=linear_model.Lasso(alpha=alpha3)
##    reg_q1=linear_model.Lasso(alpha=alpha1)
##    reg_q2=linear_model.Lasso(alpha=alpha2)
##    reg_q3=linear_model.Lasso(alpha=alpha3)
##    reg1.fit(lib1,pes[:,0])
##    reg2.fit(lib2,pes[:,1])
##    reg3.fit(lib3,pes[:,2])
##    coeff1=reg1.coef_
##    coeff2=reg2.coef_
##    coeff3=reg3.coef_
##    predict1=reg1.predict(lib1)
##    predict2=reg2.predict(lib2)
##    predict3=reg3.predict(lib3)
##    lib_q1=lib1.copy()
##    lib_q2=lib2.copy()
##    lib_q3=lib3.copy()
##    libname_q1=libname.copy()
##    libname_q2=libname.copy()
##    libname_q3=libname.copy()
##    ## choose the iteration training data to be denoisy or the ground truth
##    p_q1=pes[:,0]
##    p_q2=pes[:,1]
##    p_q3=pes[:,2]
####    p_q1=p[:,0]
####    p_q2=p[:,1]
####    p_q3=p[:,2]
##    coeff_q1=coeff1
##    coeff_q2=coeff2
##    coeff_q3=coeff3
##
##    """
##    sparse filter iterations
##    features with low correlation are deleted. The size
##    of feature library and training data is resized after
##    every iterations.
##    "fr" is the pass rate
##    "N" number of features are going to be deleted which
##        has non-zero regressed coefficient,features with
##        zero coefficient is deleted as well.
##    """
##    for iter in range(turns):
##        print(' Now processing iterration:%d'%(iter+1))
##        length_q1=coeff_q1.shape[0]
##        length_q2=coeff_q2.shape[0]
##        length_q3=coeff_q3.shape[0]
##        fr1=0.92
##        fr2=0.89
##        fr3=0.89
##        N_q1=int(np.floor((1-fr1)*length_q1))
##        N_q2=int(np.floor((1-fr2)*length_q2))
##        N_q3=int(np.floor((1-fr3)*length_q3))
##        print(length_q1,length_q2,length_q3)
##        print(N_q1,N_q2,N_q3)
##        coeffq1_sort=abs(coeff_q1)
##        coeffq2_sort=abs(coeff_q2)
##        coeffq3_sort=abs(coeff_q3)
##        coeffq1_sort.sort()
##        coeffq2_sort.sort()
##        coeffq3_sort.sort()
##        filter_index=[]
##
##        while coeffq1_sort[N_q1-1]==0:
##            N_q1=N_q1+1
##        filter_index.append(coeffq1_sort[N_q1])
##        while coeffq2_sort[N_q2-1]==0:
##            N_q2=N_q2+1
##        filter_index.append(coeffq2_sort[N_q2])
##        while coeffq3_sort[N_q3-1]==0:
##            N_q3=N_q3+1
##        filter_index.append(coeffq3_sort[N_q3])
##        q1_delete=[]
##        q2_delete=[]
##        q3_delete=[]
##        for i in range(length_q1):
##            if abs(coeff_q1[i])<filter_index[0]:
##                q1_delete.append(i)
##        for i in range(length_q2):
##            if abs(coeff_q2[i])<filter_index[1]:
##                q2_delete.append(i)
##        for i in range(length_q3):
##            if abs(coeff_q3[i])<filter_index[2]:
##                q3_delete.append(i)
##        if length_q1>11:                  
##            lib_q1=np.delete(lib_q1,q1_delete,axis=1)
##            libname_q1=np.delete(libname_q1,q1_delete,axis=0)
##        if length_q2>11:
##            lib_q2=np.delete(lib_q2,q2_delete,axis=1)
##            libname_q2=np.delete(libname_q2,q2_delete,axis=0)
##        if length_q3>8:
##            lib_q3=np.delete(lib_q3,q3_delete,axis=1)
##            libname_q3=np.delete(libname_q3,q3_delete,axis=0)    
##        reg_q1.fit(lib_q1,p_q1)
##        reg_q2.fit(lib_q2,p_q2)
##        reg_q3.fit(lib_q3,p_q3)
##        coeff_q1=reg_q1.coef_
##        coeff_q2=reg_q2.coef_
##        coeff_q3=reg_q3.coef_
##        predict_q1=reg_q1.predict(lib_q1)
##        predict_q2=reg_q2.predict(lib_q2)
##        predict_q3=reg_q3.predict(lib_q3)
##        mape1=mape(abs(p[:,0]),abs(predict1))
##        mape2=mape(abs(p[:,1]),abs(predict2))
##        mape3=mape(abs(p[:,2]),abs(predict3))
##        mse1=mse(p[:,0],predict1)
##        mse2=mse(p[:,1],predict2)
##        mse3=mse(p[:,2],predict3)
##        # MAPE is Mean Absolute Percentage Error
##        # MSE is Mean Squared Error
##        print('MAPE at τ1 after iteration(%d) is:%.10f'%(iter+1,mape1))
##        print('MSE at τ1 after iteration(%d) is:%.10f'%(iter+1,mse1))
##        print('MAPE at τ2 after iteration(%d) is:%.10f'%(iter+1,mape2))
##        print('MSE at τ2 after iteration(%d) is:%.10f'%(iter+1,mse2))
##        print('MAPE at τ3 after iteration(%d) is:%.10f'%(iter+1,mape3))
##        print('MSE at τ3 after iteration(%d) is:%.10f'%(iter+1,mse3))
##    # looking for the best alpha for each joint
##    if al in range(iter_alpha): 
##        mapelist[:,al]=mape1,mape2,mape3
##        alphalist[:,al]=alpha1,alpha2,alpha3
##        bestalpha_index=np.argmin(mapelist[0,:]),np.argmin(mapelist[1,:]),np.argmin(mapelist[2,:])
##        bestalpha=alphalist[0,bestalpha_index[0]],alphalist[1,bestalpha_index[1]],alphalist[2,bestalpha_index[2]]
##        
##    else:
##        print(coeff_q1)
##        print(libname_q1)
##        print(coeff_q2)
##        print(libname_q2)
##        print(coeff_q3)
##        print(libname_q3)
##        t1_legend.append("alpha1=%s"%(bestalpha[0]))
##        t2_legend.append("alpha2=%s"%(bestalpha[1]))
##        t3_legend.append("alpha3=%s"%(bestalpha[2]))
##        
##        # N*3 matrix denote the final coefficients for three joints at the last round
##        Xi=np.zeros([len(libname),3])
##        for i in range(len(libname_q1)):
##            for j in range(len(libname)):
##                if libname_q1[i]==libname[j]:
##                    Xi[j,0]=coeff_q1[i]
##        for i in range(len(libname_q2)):
##            for j in range(len(libname)):
##                if libname_q2[i]==libname[j]:
##                    Xi[j,1]=coeff_q2[i]
##        for i in range(len(libname_q3)):
##            for j in range(len(libname)):
##                if libname_q3[i]==libname[j]:
##                    Xi[j,2]=coeff_q1[3]
### ======================================================================================================================
##
##"""
##figures plot
##(1)measurement vs denoisy
##(2)τ1 vs sample
##(3)τ2 vs sample
##(4)τ3 vs sample
##(5)error vs alpha1
##(6)error vs alpha2
##(7)error vs alpha3
##"""
##plt.figure(2)
####plt.plot(p[:,0],'k-',label='real')
##plt.plot(pes[:,0],linewidth=0.5,color='r',label='denoisy')
##plt.plot(predict_q1,'k--',linewidth=1,label='predict_1')
##plt.legend(t1_legend)
##plt.xlabel("sample indice")
##plt.ylabel("τ1(N/m)")
##plt.grid()
##
##plt.figure(3)
####plt.plot(p[:,1],'k-',label='real')
##plt.plot(pes[:,1],linewidth=0.5,color='r',label='denoisy')
##plt.plot(predict_q2,'k--',linewidth=1,label='predict_2')
##plt.legend(t2_legend)
##plt.xlabel("sample indice")
##plt.ylabel("τ2(N/m)")
##plt.grid()
##
##plt.figure(4)
####plt.plot(p[:,2],'k-',label='real')
##plt.plot(pes[:,2],linewidth=0.5,color='r',label='denoisy')
##plt.plot(predict_q3,'k--',linewidth=1,label='predict_3')
##plt.legend(t3_legend)
##plt.xlabel("sample indice")
##plt.ylabel("τ3(N/m)")
##plt.grid()

##plt.figure(5)
##plt.plot(alphalist[0,:],mapelist[0,:])
##plt.xlabel("alpha for τ1")
##plt.ylabel("MAPE scores")
##plt.grid()
##
##plt.figure(6)
##plt.plot(alphalist[1,:],mapelist[1,:])
##plt.xlabel("alpha for τ2")
##plt.ylabel("MAPE scores")
##plt.grid()
##
##plt.figure(7)
##plt.plot(alphalist[2,:],mapelist[2,:])
##plt.xlabel("alpha for τ3")
##plt.ylabel("MAPE scores")
##plt.grid()

plt.show()



