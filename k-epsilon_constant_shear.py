#Same as model_10.py, but this one has interactive features for the user 
#input, and saves the figures
#TO BE USED FOR VALIDATION PURPOSES!


import numpy as np
import sys
import math
import matplotlib.pyplot as plt


dU_dy = 0.5

tf = 10/dU_dy
tn = tf*10.0
print tf,tn
raw_input()
#tdelta = tf/tn

Sijs = []
Ss = []
vts = []
Pks = []
ks = []
es = []
ts = []
grads = []
kolds = []
eolds = []
knews = []
enews = []
kdiffs = []
ediffs = []
iters  = []
a12s = []
a120 = 0.0
a12s.append(a120)
a11s = []
a110 = 0.0
a11s.append(a110)
iters.append(0)
################
#First order implicit scheme for the calculation of k an epsilon
################
tdelta = tf/tn

################
#Initial variables
k0 = 10.0
ks.append(k0)
eta = 2.5
residual = 0.000001
nmax = 1000000
cmu=0.09
ce1 = 1.44
ce2 = 1.92


################
#Time step 0
t0 = 0
ts.append(t0)
dU_dy0 = dU_dy#oscillating strain
grads.append(dU_dy0)

Sij0 = 0.5*( dU_dy0)
Sijs.append(Sij0)
S0 = np.sqrt(2*Sij0*Sij0)
Ss.append(S0)

e0 = (S0*k0)/eta
en = e0
#print e0
es.append(e0)

vt0 = 2*cmu*(k0**2/e0)
vts.append(vt0)
Pk0 = 2*cmu*((k0**2)/e0)*(dU_dy0**2)
Pks.append(Pk0)


kn = k0
kolds.append(kn)
en = e0
eolds.append(en)

knew = k0
knews.append(knew)
enew = e0
enews.append(enew)

kdiff = knew - k0
kdiffs.append(kdiff)
ediff = enew- e0
ediffs.append(ediff)

###################################################################
###Iterations for implicit scheme
for t in np.linspace(tdelta, tf, int(tn)):
        ts.append(t)
        grad = dU_dy
        Sij = 0.5*(grad)  #currently modified for a shear oscillating flow!
        S = np.sqrt(2*Sij*Sij)
        grads.append(grad)
        Sijs.append(Sij)
        Ss.append(S)
        #print t, grad, Sij, S

kn = k0
en = e0
k_n1_old = k0
e_n1_old = en
Pkn = 0.0
for n in range(1,len(ts)):
        for m in range(0,nmax):
                print n
                Pk_n1 = 2*cmu*((k_n1_old**2)/e_n1_old)*((Ss[n])**2)
                k_n1_new = kn + (tdelta)*(Pk_n1 - e_n1_old)
                e_n1_new = en + (tdelta)*(((Pk_n1*ce1)-(e_n1_old*ce2))*(e_n1_old/k_n1_old)) 
                
                kdiff = k_n1_new - k_n1_old
                ediff = e_n1_new - e_n1_old

                if abs(kdiff) > residual:
#                       print 'Solution is not converged!!!!!iteration no',m,'residual',abs(kdiff)
                        m += 1
                        k_n1_old = k_n1_new
                        e_n1_old = e_n1_new
                elif abs(kdiff) <= residual:
#                       print 'Solution at time step',ts[n],'is converged! Only took',m, 'iterations...'
                        iters.append(m)
                        ks.append(k_n1_new)
                        es.append(e_n1_new)
                        Pks.append(Pk_n1)
#                       print Pk_n1,k_n1_new, e_n1_new
                        vt = 2*cmu*((k_n1_new**2)/e_n1_new)
                        vts.append(vt)
                        a11 = (2.0/3.0)*k_n1_new
                        a12 = - vt*Sijs[n]/k_n1_new
                        a11s.append(a11)
                        a12s.append(a12)
                        break
        Pkn = Pk_n1     
        kn = k_n1_new
        en = e_n1_new   
        k_n1_old = k_n1_new
        e_n1_old = e_n1_new
                
               
for n in range(0,(1+int(tn))):
        item1 = ts[n]
        item2 = grads[n]
        item3 = Sijs[n]
        item4 = Ss[n]
        item5 = vts[n]
        item6 = Pks[n]
        item7 = ks[n]
        item8 = es[n]   
        item9 = iters[n]

#######Plotting 

tss = []

for n in range(0,int(len(ts))):
        s=ts[n]
        ss = s*dU_dy
        tss.append(ss)
        #print ss

ks_scaleds = []
for n in range(0,int(len(ks))):
        k_s = ks[n]
        k_scaled = k_s/k0
        ks_scaleds.append(k_scaled)
        
fig, ((ax1, ax2, ax3, ax4),(ax5, ax6, ax7, ax8)) = plt.subplots(nrows=2, ncols=4)

plt.title('Evolution of k and epsilon for homogenous cyclic strain')


ax1.plot(ts,grads,'-', color='#dc143c',linewidth = 2.0)
ax1.grid(b=True,which='both',color='0.25',linestyle='--')
ax1.set_xlabel('t ')
ax1.set_ylabel(' S= dU/dy ', fontsize = 14.0)

ax2.plot(ts,Ss,'-', color='#c71585',linewidth = 2.0)
ax2.grid(b=True,which='both',color='0.25',linestyle='--')
ax2.set_xlabel('s ')
ax2.set_ylabel('|S|', fontsize = 14.0)

ax3.plot(ts,iters,'-', color='#ff4500',linewidth = 2.0)
ax3.grid(b=True,which='both',color='0.25',linestyle='--')
ax3.set_xlabel('t ')
ax3.set_ylabel('no of iterations', fontsize = 14.0)

ax4.plot(ts,es,'-', color='#ffd700',linewidth = 2.0)
ax4.grid(b=True,which='both',color='0.25',linestyle='--')
ax4.set_xlabel('t')
ax4.set_ylabel('e', fontsize = 14.0)
        
ax5.plot(ts,Pks,'-', color='#32cd32',linewidth = 2.0)
ax5.grid(b=True,which='both',color='0.25',linestyle='--')
ax5.set_xlabel('t')
ax5.set_ylabel('Production of k', fontsize = 14.0)

ax6.plot(ts,ks,'-', color='#00ffff',linewidth = 2.0)
ax6.grid(b=True,which='both',color='0.25',linestyle='--')
ax6.set_xlabel('t ')
ax6.set_ylabel('k (turbulent kinetic energy)', fontsize = 14.0)

ax7.plot(tss,ks_scaleds,'-', color='#00008b',linewidth = 2.0)
ax7.grid(b=True,which='both',color='0.25',linestyle='--')
ax7.set_xlabel('St ')
ax7.set_ylabel('k/k0', fontsize = 14.0)

ax8.plot(tss,a11s,'-', color='#4b0082',linewidth = 2.0)
ax8.plot(tss,a12s,'--', color='#4b0082',linewidth = 2.0)
ax8.grid(b=True,which='both',color='0.25',linestyle='--')
ax8.set_xlabel('St')
ax8.set_ylabel('a12', fontsize = 14.0)

        
plt.show()


######################################################################################
#	Validation
###########################################################
#Read data from SKE 
#print '-----------------------------------------------------'
h = open('Hamlington_2008_const_shear.csv')
data3 = h.read()

#Descriptive lines
first_lines3 = data3.split(',')[:]
cc = []
k_Ham = []
St_Ham = []

for line in first_lines3:
	cc.append(line)

k_Ham = cc[0:len(first_lines3):2]
St_Ham = cc[1:(len(first_lines3)-1):2]
k_Ham.pop()


h.close()

###########################################################
#Read data from NKE 
#print '-----------------------------------------------------'
g = open('Bardina_1983_k.csv')
data2 = g.read()

#Descriptive lines
first_lines2 = data2.split(',')[:]

bb = []
k_Bar = []
St_Bar = []

for line in first_lines2:
	bb.append(line)

k_Bar = bb[0:len(first_lines2):2]
St_Bar = bb[1:(len(first_lines2)-1):2]
k_Bar.pop()


g.close()
#######################################################
#Figure plotting
fig, (ax1) = plt.subplots(nrows=1, ncols=1)


plt.title('Validation of constant shear')


ax1.plot(k_Ham,St_Ham,'-', color='#dc143c',linewidth = 2.0,label='Hamlington SKE')
ax1.plot(k_Bar,St_Bar,'o', color='black',markersize = 10.0, label='Bardina LES')
ax1.plot(tss,ks_scaleds,'-', color='#008000',linewidth = 2.0, label='My results')
ax1.grid(b=True,which='both',color='0.25',linestyle='--')
ax1.set_xlabel('St')
ax1.set_ylabel('k/k0', fontsize = 14.0)

plt.legend(loc='best')
	
plt.show()














