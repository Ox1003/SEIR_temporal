# -*- coding: utf-8 -*-
'''
SEIR_meta_distance

simulates SEIR epidemic on metapopulation network and plots arrival dates against 
effective distance values stored in matrix D
'''

import numpy as np
Metacontacts=np.load('Meta.npy')		  # contacts between classes
Metacontacts=Metacontacts + Metacontacts.T
P= Metacontacts* 1./(np.sum(Metacontacts,axis=1)) # matrix between classes
classvec= metadata[:,1]          		  # classes of the individuals
classes=np.array(['1A','1B','2A','2B','3A','3B','4A','4B','5A','5B','Teachers'])  
U= np.sum(A+B,axis=2); C=np.append(A,B,axis=2)
contacts20= np.ceil(np.sum(U)/(2*26.))      # 699 20sec-contacts per 20min window
encounters=np.sum(C!=0)/52.; avgtime=contacts20/encounters # 3.45

def transmit(s, ias,isy,home,arrival):    
    # pick random individual
    r=np.random.randint(242, size=np.ceil(encounters))
    for i in r:
        c1=np.argwhere(classes==classvec[i])[0] # class of randomly picked 1st indiv.
        # pick class of contact partner according to inter-class contact probabilty
        c2= np.argmin(np.abs(np.cumsum(P[:,c1]) - np.random.random()))
        nopupilsc2= sum(classvec==classes[c2])  # pupils in class of 2nd indiv.
        jj=np.random.randint(nopupilsc2)
        j=np.argwhere(classvec==classes[c2])[jj][0] # 2nd indiv. pick
        while j==i: 		# no contact with one self, pick new
            jj=np.random.randint(nopupilsc2); 
	    j=np.argwhere(classvec==classes[c2])[jj][0]     
        # define conditions
        infectious1= (ias[i]!=0 or isy[i]!=0) and home[i]<=0 and homecl[i]<=0 
        sus1= (s[i]==1 and homecl[i]<=0)
	sus2= (s[j]==1 and homecl[j]<=0)
        infectious2= (ias[j]!=0 or isy[j]!=0) and home[j]<=0 and homecl[j]<=0 
        if (infectious1 and sus2):                
            rtemp=np.random.random()
            if isy[i]!=0:   bina=rtemp > (1-beta)**avgtime       
            else: 	    bina=rtemp > (1-beta/2.)**avgtime
            if bina==1:
                ex[j]=1*(np.random.normal()*mu*0.1 + mu)
                arrival[j]=tt; s[j]=0   # now at exposed 
        if (infectious2 and sus1):
            rtemp=np.random.random()
            if isy[j]!=0:   bina=rtemp > (1-beta)**avgtime
            else:           bina=rtemp > (1-beta/2.)**avgtime	
            if bina==1:
                ex[i]=1*(np.random.normal()*mu*0.1 + mu)
                arrival[i]=tt; s[i]=0   # now at exposed          
    
def eod(grades,tt,isy,home,rec):            
    new=np.setdiff1d(np.nonzero(isy)[0],np.nonzero(home)[0])
    home[new]=isy[new]
    
def timepass(grades,gradesPrev,tt,ias,rec,ex,isy,home,s):
    # community infection of susceptibles
    comminf(s,ex)
    
    # those without symptoms
    ias[np.nonzero(ias)[0]]-=1.00000001
    isy[np.nonzero(isy)[0]]-=1.00000001
    # those at home
    home[np.nonzero(home)[0]]-= 1.00000001
    homecl[np.nonzero(homecl)[0]]-= 1.00000001
    
    rec[ias<0]=1
    rec[isy<0]=1
    if any(isy<0):
        for i in np.where(isy<0)[0]:
            try:
                grade=metadata[i, 1].astype('S1').astype('i1')  -1  
            except:
                grade=5 # teacher
            gradesPrev[grade,tt:slots,batchno]=gradesPrev[grade,tt,batchno]-1        
    if any(ias<0): 
        for i in np.where(ias<0)[0]:
            try:
                grade=metadata[i, 1].astype('S1').astype('i1')  -1
            except:
                grade=5     # teacher
            gradesPrev[grade,tt:slots,batchno]=gradesPrev[grade,tt,batchno]-1
    ias[ias<0]=0
    isy[isy<0]=0    

    # those exposed
    ex[np.nonzero(ex)[0]]-=1.00000001
    for i in np.argwhere(ex<0):
        
        ran=np.random.random()
        ias[i]=(ran<pa)* (np.random.normal()*gamma*0.1 + gamma)
        isy[i]=(ran>pa)*(np.random.normal()*gamma*0.1 + gamma)
        ex[i]=0
        # note: matrix A and ias, isy are already sorted 
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
        for i in range(243):
            try:
                grade=metadata[i, 1].astype('S1').astype('i1') -1
            except:
                pass  ## teachers still go to school
            if grade==grademax:
                #s[i]=0; ias[i]=0; isy[i]=0; ex[i]=0
                homecl[i]=clPeriod  

    
def comminf(s,ex):
        # infection from community
    new= np.random.rand(len(np.argwhere(s==1))) < beta_spon
    ex[np.argwhere(s==1)[new]]=1*(np.random.normal()*gamma*0.1 + gamma)    
    s[np.argwhere(s==1)[new]]=0

def svinfo(t,s_time, ias_time, isy_time, rec_time,ex_time, s,ias,isy,rec,ex):
	# svinfo saves the current states and gives t=t+1
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
    maxbatch=20
    maxdays=50;  noplot=0
    threshhold=3; clPeriod=24*3*4;
    slots=24*3*maxdays +25*3*np.ceil(maxdays/7.)
    epd=[]
    
    gradesPrev=np.zeros((6,slots,maxbatch),dtype=np.int)
    for batchno in range(maxbatch):
        s_time= np.zeros(slots,dtype=np.int)
        ias_time= np.zeros(slots,dtype=np.int)
        isy_time= np.zeros(slots,dtype=np.int)
        rec_time= np.zeros(slots,dtype=np.int)
        ex_time= np.zeros(slots,dtype=np.int)
        home=np.zeros((no),dtype=np.float)
        homecl=np.zeros((no),dtype=np.float)
        s=np.ones((no),dtype=np.int)
        ias=np.zeros((no),dtype=np.float)
        isy=np.zeros((no),dtype=np.float)
        ex=np.zeros((no),dtype=np.float)
        rec=np.zeros((no),dtype=np.int)
        previousgrade= np.zeros(5,dtype=np.int)
        classes=np.array(['1A','1B','2A','2B','3A','3B','4A','4B','5A','5B','Teachers']) 
        
        arrival=np.zeros(242)
        tt=0;
        grades=np.zeros((6,slots),dtype=np.int)        
        gamma=24*3*4.
        mu=24*3*2.
        beta=20*3.5* 10**-4 ; beta_spon=20*2.8* 10**-9
        pa=0.3 
        dates=[]
        tinder=np.random.randint(no, size=1)
        isy[tinder]=1*(np.random.normal()*gamma*0.1 + gamma)
        s[tinder]=0
        print tinder
        # initial
        try:
            grade=metadata[tinder, 1].astype('S1').astype('i1')  -1  
        except:
            grade=5 # teacher       
        gradesPrev[grade,tt:slots,batchno]=gradesPrev[grade,tt,batchno]+1  
        while (tt<maxdays*24*3):
            if not any(gradesPrev[:,tt,batchno]) and all(ex==0):
                break
            else:
                if np.mod(tt+24*3*2,24*3*7)==0:    #weekend 
                    for t in range (72+72):
                        timepass(grades,gradesPrev,tt,ias,rec,ex,isy,home,s)
                     #    svinfo(tt,s_time, ias_time, isy_time, rec_time,ex_time, s,ias,isy,rec,ex)
                        tt+=1
                    eod(grades,tt,isy,home,rec)
                for t in day1:
                    transmit(s,ias,isy,home,arrival)
                    timepass(grades,gradesPrev,tt,ias,rec,ex,isy,home,s)
                  #  targetedClosure(tt,gradesPrev,ias,isy,home,homecl,s,dates)
                  #   svinfo(tt,s_time, ias_time, isy_time, rec_time,ex_time, s,ias,isy,rec,ex)
                    tt+=1                     
                for t in range (72-25):
                    timepass(grades,gradesPrev,tt,ias,rec,ex,isy,home,s)
                    svinfo(tt,s_time, ias_time, isy_time, rec_time,ex_time, s,ias,isy,rec,ex)
                    tt+=1
                eod(grades,tt,isy,home,rec)
                if np.mod(tt+24*3*4,24*3*7)==0:    #wednesday
                    for t in range (72):
                        timepass(grades,gradesPrev,tt,ias,rec,ex,isy,home,s)
                    #     svinfo(tt,s_time, ias_time, isy_time, rec_time,ex_time, s,ias,isy,rec,ex)
                        tt+=1
                    eod(grades,tt,isy,home,rec)
                if np.mod(tt+24*3*2,24*3*7)==0:    #weekend 
                    for t in range (72+72):
                        timepass(grades,gradesPrev,tt,ias,rec,ex,isy,home,s)
                   #      svinfo(tt,s_time, ias_time, isy_time, rec_time,ex_time, s,ias,isy,rec,ex)
                        tt+=1 
                    eod(grades,tt,isy,home,rec)
                for t in day2:
                    transmit(s,ias,isy,home,arrival)
                    timepass(grades,gradesPrev,tt,ias,rec,ex,isy,home,s)
                 #   targetedClosure(tt,gradesPrev,ias,isy,home,homecl,s,dates)
                    svinfo(tt,s_time, ias_time, isy_time, rec_time,ex_time, s,ias,isy,rec,ex)
                    tt+=1       
                for t in range (72-25):
                    timepass(grades,gradesPrev,tt,ias,rec,ex,isy,home,s)
                    svinfo(tt,s_time, ias_time, isy_time, rec_time,ex_time, s,ias,isy,rec,ex)
                    tt+=1
                eod(grades,tt,isy,home,rec)
                if np.mod(tt+24*3*4,24*3*7)==0:    #wednesday
                    for t in range (72):
                        timepass(grades,gradesPrev,tt,ias,rec,ex,isy,home,s)
                        svinfo(tt,s_time, ias_time, isy_time, rec_time,ex_time, s,ias,isy,rec,ex)
                        tt+=1
                    eod(grades,tt,isy,home,rec)

        #############################################################
        numInfec=np.sum(gradesPrev[:,:,batchno], axis=0)    
        try:   # sometimes error occurs when data too small for regression
            D=np.load('D.npy')
            plt.close('all')
            classes=np.array(['1A','1B','2A','2B','3A','3B','4A','4B','5A','5B','Te']) 
            tinderstring=metadata[tinder, 1].astype('S2')
            tinderclass= np.argwhere(classes==tinderstring)[0][0]
            arrivalclass=np.zeros((11,243))
            d4=np.zeros(242)
            for i in np.nonzero(arrival)[0]:
                a=metadata[i, 1].astype('S2') 
                arrivalclass[np.argwhere(classes==a)[0],i]=  arrival[i]
                d4[i]=D[tinderclass,np.argwhere(classes==a)[0][0]]
            arrivalclass=np.max(arrivalclass, axis=0)
            fig=plt.figure()
            plt.rcParams.update({'font.size': 14})
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif')
            y=arrivalclass.T/72.
            plt.subplot(1,1,1)  
            plt.title(r'$\mbox{Effective distances vs date of first infection$')
            x=d4 
            plt.scatter(x[x!=0], y[x!=0],marker='x',c='darkblue')    
            fit =np.polyfit(x[x!=0], y[x!=0],1)
            slope, _, r_value, _, std_err = stats.linregress(x[x!=0], y[x!=0])
            z = np.poly1d(fit)
            plt.plot(x, z(x),'k-', linewidth=1)
            plt.legend([r'$\mbox{linear fit, } R = %s$'%(np.round(r_value,2)),'Data'],fontsize=12)
            plt.ylabel('Exposure date')
            plt.xlabel('Effective class distance')
            plt.xlim(xmin=0)
            plt.ylim(ymin=0)
            plt.show()
            fig.savefig('effDistMETADt_bta3.5_%s.png'%(tinder), dpi=200) 
        except:
            pass
        
