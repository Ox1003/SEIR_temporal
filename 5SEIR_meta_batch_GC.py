# -*- coding: utf-8 -*-
'''
SEIR_meta_batch_GC

simulates SEIR epidemic on the metapopulation network. Can specify batch size.
Here, targeted grade closure (GC) is implemented.

@author: 
EnCt28648cf14f7c30972de440f01c53b7544e83c9ee78648cf14f7c30972de440f01puAufN5B6QB
rpWuQm1ZPfSb3bcvooZYGIwEmS

Decrypt it at https://encipher.it
'''
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

A=np.load('day1.npy')
B=np.load('day2.npy')
no=np.load('numbers.npy') +1
day1=np.load('times1.npy')
day2=np.load('times2.npy')
metadata=np.loadtxt('metadata_primaryschool.txt',delimiter='\t',dtype='str')
Meta=np.load('Meta.npy')
metadata=np.loadtxt('metadata_primaryschool.txt',delimiter='\t',dtype='S16')
classes=np.array(['1A','1B','2A','2B','3A','3B','4A','4B','5A','5B','Teachers'])  
Metacontacts=np.load('Meta.npy')
Metacontacts=Metacontacts + Metacontacts.T
P= Metacontacts* 1./(np.sum(Metacontacts,axis=1))      # koennte auch dynamische P matrix nehmen ..
classvec= metadata[:,1]
U= np.sum(A+B,axis=2)
C=np.append(A,B,axis=2)
#159 contacts per 20min window during school hours
contacts20= np.ceil(np.sum(U)/(2*26.))   #2400 20sec-contacts per 20min window during school hours
encounters=np.sum(C!=0)/52.
avgtime=contacts20/encounters
def transmit(s, ias,isy,home,arrival):
    
    r=np.random.randint(242, size=np.ceil(encounters))
   # for i in r[(isy+ias)[r]!=0]:
    for i in r:
        c1=np.argwhere(classes==classvec[i])[0][0]
        c2= np.argmin(np.abs(np.cumsum(P[:,c1]) - np.random.random()))
        nopupilsc2= sum(classvec==classes[c2])
        jj=np.random.randint(nopupilsc2)
        j=np.argwhere(classvec==classes[c2])[jj][0]
        while j==i:
            jj=np.random.randint(nopupilsc2)
            j=np.argwhere(classvec==classes[c2])[jj][0]
        infectious1= (ias[i]!=0 or isy[i]!=0)and home[i]<=0 and homecl[i]<=0 
        sus1= s[i]==1 and homecl[i]<=0 
        infectious2=(ias[j]!=0 or isy[j]!=0)and home[j]<=0 and homecl[j]<=0 
        sus2= s[j]==1 and homecl[j]<=0
        if (infectious1 and sus2):                
            rtemp=np.random.random()
           # avgtime=C[np.nonzero(C)][np.random.randint(len(C))]
            if isy[i]!=0:
                bina=rtemp > (1-beta)**avgtime
            else:
                bina=rtemp > (1-beta/2.)**avgtime
            if bina==1:
               # print 'pupil w/ ID', j, 'got exposed'
                ex[j]=1*(np.random.normal()*mu*0.1 + mu)
                arrival[j]=tt
                s[j]=0   #now at exposed 
        if (infectious2 and sus1):
            rtemp=np.random.random()
            #avgtime=C[np.nonzero(C)][np.random.randint(len(C))]
            if isy[j]!=0:
                bina=rtemp > (1-beta)**avgtime
            else:
                bina=rtemp > (1-beta/2.)**avgtime
            if bina==1:
               # print 'pupil w/ ID', j, 'got exposed'
                ex[i]=1*(np.random.normal()*mu*0.1 + mu)
                arrival[i]=tt
                s[i]=0   #now at exposed
                       
    
def eod(grades,tt,isy,home,rec):            
    new=np.setdiff1d(np.nonzero(isy)[0],np.nonzero(home)[0])
    home[new]=isy[new]
    
def timepass(grades,gradesPrev,tt,ias,rec,ex,isy,home,s):
    #spontaneous infection of susceptibles
    spontexpo(s,ex)
    
    #those without symptoms
    ias[np.nonzero(ias)[0]]-=1.00000001
    isy[np.nonzero(isy)[0]]-=1.00000001
    #those at home
    home[np.nonzero(home)[0]]-= 1.00000001
    homecl[np.nonzero(homecl)[0]]-= 1.00000001
    
    rec[ias<0]=1
    rec[isy<0]=1
    if any(isy<0):
        for i in np.where(isy<0)[0]:
            try:
                grade=metadata[i, 1].astype('S1').astype('i1')  -1  
            except:
                grade=5 #teacher
            gradesPrev[grade,tt:slots,batchno]=gradesPrev[grade,tt,batchno]-1        
    if any(ias<0): 
        for i in np.where(ias<0)[0]:
            try:
                grade=metadata[i, 1].astype('S1').astype('i1')  -1
            except:
                grade=5     #teacher
            gradesPrev[grade,tt:slots,batchno]=gradesPrev[grade,tt,batchno]-1
    ias[ias<0]=0
    isy[isy<0]=0    
    
    

    #those exposed
    ex[np.nonzero(ex)[0]]-=1.00000001
    for i in np.argwhere(ex<0):
        
        ran=np.random.random()
        ias[i]=(ran<pa)* (np.random.normal()*gamma*0.1 + gamma)
        isy[i]=(ran>pa)*(np.random.normal()*gamma*0.1 + gamma)
        ex[i]=0
        #note: my matrix A and ias, isy are already sorted 
        try:
            grade=metadata[i, 1].astype('S1').astype('i1')[0]  -1
            grades[grade,tt:slots]=grades[grade,tt]+1
            gradesPrev[grade,tt:slots,batchno]=gradesPrev[grade,tt,batchno]+1
        except:
            gradesPrev[5,tt:slots,batchno]=gradesPrev[5,tt,batchno]+1   
    
def targetedClosure(tt,gradesPrev,ias,isy,home,homecl,s,dates):
    grademax=np.argmax(gradesPrev[:,tt,batchno]) 
    if gradesPrev[grademax,tt,batchno] > threshhold and previousgrade[grademax]==0:
        dates.append(tt/72.0)
        previousgrade[grademax]=1
        #gradesPrev[grademax,tt:slots]=05
        for i in range(243):
            try:
                grade=metadata[i, 1].astype('S1').astype('i1') -1
            except:
                pass  ## teachers still go to school
            if grade==grademax:
                #s[i]=0; ias[i]=0; isy[i]=0; ex[i]=0
                homecl[i]=clPeriod  

    
def spontexpo(s,ex):
	#spontinfec accounts for the spontaenous exposure
    new= np.random.rand(len(np.argwhere(s==1))) < beta_spon
  #  if max(new)!=0:
  #      print 'spont infection of', np.argwhere(s==1)[new]
    ex[np.argwhere(s==1)[new]]=1*(np.random.normal()*gamma*0.1 + gamma)    
    s[np.argwhere(s==1)[new]]=0

def svinfo(t,s_time, ias_time, isy_time, rec_time,ex_time, s,ias,isy,rec,ex):
	## svinfo saves the current states and gives t=t+1
    s_time[t]=(sum(s!=0))
    ias_time[t]=(sum(ias!=0))
    isy_time[t]=(sum(isy!=0))
    rec_time[t]=(sum(rec!=0))
    ex_time[t]=(sum(ex!=0))
        
    
if __name__ == "__main__": 
    fig=plt.figure()
    plt.rcParams.update({'font.size': 12})
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    maxbatch=100
    maxdays=40; tt=0; 
    threshhold=3; clPeriod=24*3*4.;
    previousgrade= np.zeros(5,dtype=np.int)
    slots=24*3*maxdays +np.ceil(maxdays/7) *24 + 1
    noplot=0; epd=[]
    gradesPrev=np.zeros((6,slots,maxbatch),dtype=np.int)
    for batchno in range(maxbatch):    
        print batchno
        home=np.zeros((no),dtype=np.float)
        homecl=np.zeros((no),dtype=np.float)
        s=np.ones((no),dtype=np.int)
        ias=np.zeros((no),dtype=np.float)
        isy=np.zeros((no),dtype=np.float)
        ex=np.zeros((no),dtype=np.float)
        rec=np.zeros((no),dtype=np.int)
        tt=0
#        s_time= np.zeros(slots,dtype=np.int)
#        ias_time= np.zeros(slots,dtype=np.int)
#        isy_time= np.zeros(slots,dtype=np.int)
#        rec_time= np.zeros(slots,dtype=np.int)
#        ex_time= np.zeros(slots,dtype=np.int)
        grades=np.zeros((6,slots),dtype=np.int)        
        classPrev=np.zeros((11,slots),dtype=np.int)
        gamma=24*3*4.
        mu=24*3*2.
        beta=20*3.5* 10**-4 ; beta_spon=20*2.8* 10**-9
        pa=0.3 
        dates=[]
        arrival=np.zeros((243,1))
        tinder=np.random.randint(no, size=1)
        isy[tinder]=1*(np.random.normal(len(tinder))*gamma*0.1 + gamma)
        s[tinder]=0
        print tinder
        #initial
        try:
            grade=metadata[tinder, 1].astype('S1').astype('i1')  -1  
        except:
            grade=5 #teacher 
           # tinderstring=metadata[tinder, 1].astype('S2')
           # tinderclass= np.argwhere(classes==tinderstring)[0][0]
        grades[grade,tt:slots]=grades[grade,tt]+1
        gradesPrev[grade,tt:slots,batchno]=gradesPrev[grade,tt,batchno]+1   
        while tt<maxdays*24*3:
            if not any(gradesPrev[:,tt,batchno]) and all(ex==0):
                break
            else:
                if np.mod(tt+24*3*2,24*3*7)==0:    #weekend 
                    for t in range (72+72):
                        timepass(grades,gradesPrev,tt,ias,rec,ex,isy,home,s)
                    #    svinfo(tt,s_time, ias_time, isy_time, rec_time,ex_time, s,ias,isy,rec,ex)
                        tt+=1 
                    eod(grades, tt,isy,home,rec)
                for t in day1:
                    transmit(s, ias,isy,home,arrival)
                    timepass(grades,gradesPrev,tt,ias,rec,ex,isy,home,s)
                    targetedClosure(tt,gradesPrev,ias,isy,home,homecl,s,dates)
                 #   svinfo(tt,s_time, ias_time, isy_time, rec_time,ex_time, s,ias,isy,rec,ex)
                    tt+=1          
                
                for t in range (72-25):
                    timepass(grades,gradesPrev,tt,ias,rec,ex,isy,home,s)
                 #   svinfo(tt,s_time, ias_time, isy_time, rec_time,ex_time, s,ias,isy,rec,ex)
                    tt+=1
                eod(grades,tt,isy,home,rec)
                if np.mod(tt+24*3*2,24*3*7)==0:    #weekend 
                    for t in range (72+72):
                        timepass(grades,gradesPrev,tt,ias,rec,ex,isy,home,s)
                #        svinfo(tt,s_time, ias_time, isy_time, rec_time,ex_time, s,ias,isy,rec,ex)
                        tt+=1 
                    eod(grades,tt,isy,home,rec)
                for t in day2:
                    transmit(s, ias,isy,home,arrival)
                    timepass(grades,gradesPrev,tt,ias,rec,ex,isy,home,s)
                    targetedClosure(tt,gradesPrev,ias,isy,home,homecl,s,dates)
                #    svinfo(tt,s_time, ias_time, isy_time, rec_time,ex_time, s,ias,isy,rec,ex)
                    tt+=1       
                for t in range (72-25):
                    timepass(grades,gradesPrev,tt,ias,rec,ex,isy,home,s)
                #    svinfo(tt,s_time, ias_time, isy_time, rec_time,ex_time, s,ias,isy,rec,ex)
                    tt+=1
                eod(grades,tt,isy,home,rec)
    #############################################################  
        print np.max(np.sum(gradesPrev[:,:,batchno], axis=0))
        if np.max(np.sum(gradesPrev[:,:,batchno], axis=0))<23   :
            #gradesPrev[:,:,batchno]=0
            print 'less than 10 percent', batchno
            noplot+=1
        else:
            epd.append(batchno)
            print 'bigger 10 perc', batchno

   # plt.plot(np.arange(slots)/72., np.sum(gradesPrev, axis=0) ,linewidth=0.5)  
    infecs=np.sum(gradesPrev, axis=0)
    plt.plot(np.arange(slots)/72., np.median(infecs.T[epd],axis=0) ,linewidth=1, color='r') 
    plt.plot(np.arange(slots)/72., np.median(infecs.T,axis=0) ,linewidth=1, color='k') 
    plt.plot(np.arange(slots)/72., infecs,linewidth=0.3,alpha=0.2)
    
    plt.title(r'$\mbox{Prevalence of SEIR disease; %s simulations, %s with AR\textgreater0.1}$'%(maxbatch, maxbatch-noplot))
    plt.xlabel('time (days)')
    plt.ylabel('Number infected')
    plt.legend([r'$\mbox{median (only AR\textgreater0.1})$', 'median'])
    
    fig.savefig('%sMETAgc_median_Infecs_bta3.5_gam4_th3_clP4avg.png'%maxbatch, bbox_inches='tight',dpi=500)
    plt.show()
