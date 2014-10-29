#Plot Data

import numpy as np
import matplotlib.pyplot as plt
import sys


validation = 0
#######################################################################################    
def data_extraction_keps(validation, casetype):

    f = open('Output_keps.dat')
    data = f.read()
    lines = data.split('\n')
    lines.pop()
    lines.pop(0)

    ns = []
    xs = []
    ts = []
    ss = []
    pks = []
    ks = []
    es = []

    for item in lines:
        #print 'item is:', item
        n = item[0:15]
        ns.append(n)
        x = item[15:18]
        xs.append(x)
        t = item[17:33]
        ts.append(t)
        s = item[34:50]
        ss.append(s)
        pk = item[50:64]
        pks.append(pk)
        k = item[64:80]
        ks.append(k)
        e = item[80:96]
        es.append(e)
        #print 'n', n
        #print 'x', x
        #print 't', t
        #print 's', s
        #print 'pk', pk
        #print 'k', k
        #print 'e', e

    #print len(ns), len(xs), len(ts), len(ss), len(pks), len(ks), len(es)

    fig, ((ax1, ax2),(ax3, ax4)) = plt.subplots(nrows=2, ncols=2)

    ax1.plot(ts,ss,'-', color='black',linewidth = 2.0)
    ax1.grid(b=True,which='both',color='0.25',linestyle='--')
    ax1.tick_params(axis = 'x', labelsize = 12.0)
    ax1.tick_params(axis = 'y', labelsize = 12.0)
    ax1.set_xlabel('t(s) ', fontsize = 12.0)
    ax1.set_ylabel(' S(t)', fontsize = 12.0)
    
    ax2.plot(ts,ks,'-', color='blue',linewidth = 2.0)
    ax2.grid(b=True,which='both',color='0.25',linestyle='--')
    ax2.tick_params(axis = 'x', labelsize = 12.0)
    ax2.tick_params(axis = 'y', labelsize = 12.0)
    ax2.set_xlabel('t (s) ', fontsize = 12.0)
    ax2.set_ylabel('ks ', fontsize = 12.0)
    
    ax3.plot(ts,pks,'-', color='red',linewidth = 2.0)
    ax3.grid(b=True,which='both',color='0.25',linestyle='--')
    ax3.tick_params(axis = 'x', labelsize = 12.0)
    ax3.tick_params(axis = 'y', labelsize = 12.0)
    ax3.set_xlabel('t(s) ', fontsize = 12.0)
    ax3.set_ylabel('Pks ', fontsize = 12.0)
    
    ax4.plot(ts,ks,'-', color='green',linewidth = 2.0)
    ax4.grid(b=True,which='both',color='0.25',linestyle='--')
    ax4.tick_params(axis = 'x', labelsize = 12.0)
    ax4.tick_params(axis = 'y', labelsize = 12.0)
    ax4.set_xlabel('t(s) ', fontsize = 12.0)
    ax4.set_ylabel('$\epsilon (s)$', fontsize = 12.0)
    
    plt.show()
    
    if int(validation) == 1:
        if int(casetype) ==  1:
            print 'Validating Results...'
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
            ax1.plot(ts,ks,'-', color='#008000',linewidth = 2.0, label='My results')
            ax1.grid(b=True,which='both',color='0.25',linestyle='--')
            ax1.set_xlabel('St')
            ax1.set_ylabel('k/k0', fontsize = 14.0)
            ax1.set_xlim([0,6.0])
            ax1.set_ylim([0,5.0])

            plt.legend(loc='best')
	
            plt.show()
        
        elif int(casetype) == 2:
            print 'No validation data exits!'
            
        elif int(casetype) == 3:
            #Hadzic & Hanjalic 
            print 'Hadzic & Hanjalic'
            
        elif int(casetype) == 4:
            ###########################################################
            #Read data from Chen2006 
            #print '-----------------------------------------------------'
            f = open('Chen_Ham2008.csv')
            data = f.read()

            #Descriptive lines
            first_lines = data.split(',')[:]
            aa = []
            Ss_Chen = []
            ts_Chen = []

            for line in first_lines:
            	aa.append(line)

            Ss_Chen = aa[0:len(first_lines):2]
            ts_Chen = aa[1:(len(first_lines)-1):2]
            Ss_Chen.pop()

            f.close()

            ###########################################################
            #Read data from NKE 
            #print '-----------------------------------------------------'
            g = open('NKE_Ham2008_1.csv')
            data2 = g.read()

            #Descriptive lines
            first_lines2 = data2.split(',')[:]

            bb = []
            Ss_NKE = []
            ts_NKE = []

            for line in first_lines2:
            	bb.append(line)

            Ss_NKE = bb[0:len(first_lines2):2]
            ts_NKE = bb[1:(len(first_lines2)-1):2]
            Ss_NKE.pop()


            g.close()

            ###########################################################
            #Read data from SKE 
            #print '-----------------------------------------------------'
            h = open('SKE_Ham2008_1.csv')
            data3 = h.read()

            #Descriptive lines
            first_lines3 = data3.split(',')[:]
            cc = []
            Ss_SKE = []
            ts_SKE = []

            for line in first_lines3:
            	cc.append(line)

            Ss_SKE = cc[0:len(first_lines3):2]
            ts_SKE = cc[1:(len(first_lines3)-1):2]
            Ss_SKE.pop()


            h.close()


          #  print len(Ss_Chen), len(ts_Chen)
          #  print len(Ss_NKE), len(ts_NKE)
          #  print len(Ss_SKE), len(ts_SKE)

          #  print Ss_NKE, ts_NKE
          #  print 
          #  print Ss_SKE, ts_SKE

            #######################################################
            #Figure plotting
            fig, (ax1) = plt.subplots(nrows=1, ncols=1)


            plt.title('Validation of time variant strain')


            ax1.plot(Ss_Chen,ts_Chen,'o', color='#dc143c',linewidth = 2.0, markersize=12.0,label='Chen at all')
            ax1.plot(Ss_NKE,ts_NKE,'*', color='#0000cd',linewidth = 2.0, label='Hamlington & Dahm NKE')
            ax1.plot(Ss_SKE,ts_SKE,'-', color='#008000',linewidth = 2.0, label='Hamlington & Dahm SKE')
            #ax1.plot(t_scaleds, a11s,'--', color='#cc66ff',linewidth = 2.0, label='My results')
            ax1.grid(b=True,which='both',color='0.25',linestyle='--')
            ax1.set_xlabel('t e0/k0')
            ax1.set_ylabel('a11', fontsize = 14.0)

            plt.legend(loc='best')
	
            plt.show()
            
        
    else:
        print 'No Validation'
        

    
    f.close()
    return
#######################################################################################    
def data_extraction_RSM(validation, casetype):

    ns = []
    xs = []
    ts = []
    ss = []
    pks = []
    ks = []
    es = []
    a11s = []
    a22s = []
    a33s = []
    a12s = []
    u1u1s = []
    u2u2s = []
    u3u3s = []
    u1u2s = []
    
    
    f = open('Output_RSM.dat')
    data = f.read()
    lines = data.split('\n')
    lines.pop()
    lines.pop(0)            

    for item in lines:
        #print 'item is:', item
        n = item[0:15]
        ns.append(n)
        x = item[15:18]
        xs.append(x)
        t = item[17:33]
        ts.append(t)
        s = item[34:50]
        ss.append(s)
        pk = item[50:64]
        pks.append(pk)
        k = item[64:80]
        ks.append(k)
        e = item[80:96]
        es.append(e)
        #print 'n', n
        #print 'x', x
        #print 't', t
        #print 's', s
        #print 'pk', pk
        #print 'k', k
        #print 'e', e

    #print len(ns), len(xs), len(ts), len(ss), len(pks), len(ks), len(es)

    f.close()
    
    g = open('Output_RSM_stresses.dat')
    data1 = g.read()
    lines1 = data1.split('\n')
    lines1.pop()
    lines1.pop(0)
    for item in lines1:
        #print item
        a11 = item[1:17]
        a11s.append(a11)
        a22 = item[18:33]
        a22s.append(a22)
        a33 = item[34:49]
        a33s.append(a33)
        a12 = item[50:65]
        a12s.append(a12)
        u1u1 = item[66:81]
        u1u1s.append(u1u1)
        u2u2 = item[82:97]
        u2u2s.append(u2u2)
        u3u3 = item[98:113]
        u3u3s.append(u3u3)
        u1u2 = item[114:129]
        u1u2s.append(u1u2)
        
        #print 'a11', a11
        #print 'a22', a22
        #print 'a33', a33
        #print 'a12', a12
        #print 'u1u1', u1u1
        #print 'u2u2', u2u2
        #print 'u3u3', u3u3
        #print 'u1u2', u1u2
        #raw_input()
        
    g.close()

    #print len(a11s), len(a22s), len(a33s), len(a12s), len(u1u1s), len(u2u2s), len(u3u3s), len(u1u2)

    fig, ((ax1, ax2),(ax3, ax4)) = plt.subplots(nrows=2, ncols=2)

    ax1.plot(ts,ss,'-', color='black',linewidth = 2.0)
    ax1.grid(b=True,which='both',color='0.25',linestyle='--')
    ax1.tick_params(axis = 'x', labelsize = 12.0)
    ax1.tick_params(axis = 'y', labelsize = 12.0)
    ax1.set_xlabel('t(s) ', fontsize = 12.0)
    ax1.set_ylabel(' S(t)', fontsize = 12.0)

    ax2.plot(ts,ks,'-', color='blue',linewidth = 2.0)
    ax2.grid(b=True,which='both',color='0.25',linestyle='--')
    ax2.tick_params(axis = 'x', labelsize = 12.0)
    ax2.tick_params(axis = 'y', labelsize = 12.0)
    ax2.set_xlabel('t (s) ', fontsize = 12.0)
    ax2.set_ylabel('ks ', fontsize = 12.0)

    ax3.plot(ts,pks,'-', color='red',linewidth = 2.0)
    ax3.grid(b=True,which='both',color='0.25',linestyle='--')
    ax3.tick_params(axis = 'x', labelsize = 12.0)
    ax3.tick_params(axis = 'y', labelsize = 12.0)
    ax3.set_xlabel('t(s) ', fontsize = 12.0)
    ax3.set_ylabel('Pks ', fontsize = 12.0)

    ax4.plot(ts,ks,'-', color='green',linewidth = 2.0)
    ax4.grid(b=True,which='both',color='0.25',linestyle='--')
    ax4.tick_params(axis = 'x', labelsize = 12.0)
    ax4.tick_params(axis = 'y', labelsize = 12.0)
    ax4.set_xlabel('t(s) ', fontsize = 12.0)
    ax4.set_ylabel('$\epsilon (s)$', fontsize = 12.0)

    plt.show()


    fig, ((ax5),(ax6)) = plt.subplots(nrows=2, ncols=1)
    
    ax5.plot(ts,a11s,'-', color='#660000',linewidth = 2.0, label='$a_{11}$')
    ax5.plot(ts,a22s,'-', color='#336600',linewidth = 2.0, label='$a_{22}$')
    ax5.plot(ts,a33s,'-', color='#bb00ff',linewidth = 2.0, label='$a_{33}$')
    ax5.plot(ts,a12s,'-', color='#ff4500',linewidth = 2.0, label='$a_{12}$')
    ax5.grid(b=True,which='both',color='0.25',linestyle='--')
    ax5.tick_params(axis = 'x', labelsize = 15.0)
    ax5.tick_params(axis = 'y', labelsize = 15.0)
    ax5.set_xlabel('t(s) ', fontsize = 15.0)
    ax5.set_ylabel('$a_{ij}$ ', fontsize = 15.0)
    ax5.legend(bbox_to_anchor=(1.005, 1), loc=2, borderaxespad=0.)

    ax6.plot(ts,u1u1s,'-', color='#660000',linewidth = 2.0, label='$u_1 u_1$')
    ax6.plot(ts,u2u2s,'-', color='#336600',linewidth = 2.0, label='$u_2 u_2$')
    ax6.plot(ts,u3u3s,'-', color='#bb00ff',linewidth = 2.0, label='$u_3 u_3$')
    ax6.plot(ts,u1u2s,'-', color='#ff4500',linewidth = 2.0, label='$u_1 u_2$')
    ax6.grid(b=True,which='both',color='0.25',linestyle='--')
    ax6.tick_params(axis = 'x', labelsize = 15.0)
    ax6.tick_params(axis = 'y', labelsize = 15.0)
    ax6.set_xlabel('t(s) ', fontsize = 15.0)
    ax6.set_ylabel('$u_i u_j$', fontsize = 15.0)
    ax6.legend(bbox_to_anchor=(1.005, 1), loc=2, borderaxespad=0.)

    plt.show()
    
    return    
#######################################################################################    

print '1 - k-epsilon; 2 - RSM, 3 - Cas'
turbmod=int(raw_input('WHAT TURBULENCE MODEL:   '))
if turbmod!=1 and turbmod!=2 and turbmod!=3:
    print 'Error'
    sys.exit()

print '0 - NO; 1 - YES'
validation = int(raw_input('VALIDATE RESULTS? '))
if validation!= 1 and validation!= 0 :
    print 'Error'
    sys.exit()
    
print '1 - A; 2 - B; 3 - C; 4 - D; 5 - E'
casetype = int(raw_input('WHAT CASE? '))
if casetype!=1 and casetype!=2 and casetype!= 3 and casetype!=4 and casetype!= 5:
    print 'Error'
    sys.exit()

if int(turbmod) == 1:
    data_extraction_keps(validation, casetype)

elif int(turbmod) == 2:
    data_extraction_RSM(validation, casetype)

elif int(turbmod) == 3:
    data_extraction_Cas(validation, casetype)

else:
    print 'I\'m going fishing!'
    
    
