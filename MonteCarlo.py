import scipy.stats as stats
import numpy as np

def LinRegress(inputs,y):
    E=np.ones((len(inputs[0]),len(inputs)+1))
    for i in range(len(inputs)):
        E[:,i]=inputs[i]
    #E[:,-1]=np.random.normal(size=E[:,-1].shape[0])
    xhat = np.dot(np.dot(np.linalg.inv(np.dot(np.transpose(E),E)),np.transpose(E)),y)
    return E,xhat

def SigTrend(data,trend,num_iter,alpha=5.0):
    print(data.shape)
    print(trend.shape)
    Y,S = data.shape
    #print(Y,I,J)
    #It,Jt = trend.shape
    #data.shape=(Y,I*J)
    #trend.shape=(It*Jt,)
    Years=np.arange(Y)
    Shuffle=np.arange(Y)
    MCPerms=np.ones((num_iter,S))*np.nan
    P=np.ones(trend.shape)*np.nan
    for i in range(num_iter):
        np.random.shuffle(Years)
        E,xhat=LinRegress([Years[Shuffle]],data)
        MCPerms[i]=xhat[0]

    for i in range(S):
        P[i]=stats.percentileofscore(MCPerms[:,i],trend[i])

    P_tf=(P>=(100-alpha/2.0)) | (P<=alpha/2.0)

    #data.shape=(I,J)
    #trend.shape=(I,J)
    #P.shape=(I,J)
    #P_tf.shape=(I,J)
    return P,P_tf

def Monte_Carlo_1D(Obs,TS,NumDraw,CompType='eq',CompVal=0,num_iter=5000):
    #For 1 dimensional data such as a time series at a single point
    #Need to modify for larger dimentional datasets
    MCPerms=np.ones((num_iter,))*np.nan
    for i in range(num_iter):
        shuffle_inds=np.arange(TS.shape[0])
        np.random.shuffle(shuffle_inds)
        if CompType=='ge':
            MCPerms[i]=(TS[shuffle_inds][:NumDraw]>=CompVal).sum()
        elif CompType=='gt':
            MCPerms[i]=(TS[shuffle_inds][:NumDraw]>CompVal).sum()
        elif CompType=='eq':
            MCPerms[i]=(TS[shuffle_inds][:NumDraw]==CompVal).sum()
        elif CompType=='lt':
            MCPerms[i]=(TS[shuffle_inds][:NumDraw]<CompVal).sum()
        elif CompType=='le':
            MCPerms[i]=(TS[shuffle_inds][:NumDraw]<=CompVal).sum()

    P=stats.percentileofscore(MCPerms,Obs)
    return P

def Monte_Carlo_2D(obs,TS,NumDraw,CompType='eq',CompVal=0,num_iter=5000):
    T,I,J=TS.shape
    MCPerms=np.ones((num_iter,I,J))
    for i in range(num_iter):
        shuffle_inds=np.arange(TS.shape[0])
        np.random.shuffle(shuffle_inds)
        if CompType=='ge':
            MCPerms[i]=(TS[shuffle_inds][:NumDraw]>=CompVal).sum(axis=0)
        elif CompType=='gt':
            MCPerms[i]=(TS[shuffle_inds][:NumDraw]>CompVal).sum(axis=0)
        elif CompType=='eq':
            print(TS[shuffle_inds])
            print(NumDraw)
            print(CompVal)
            print(i)
            print(MCPerms.shape)
            print((TS[shuffle_inds][:NumDraw]==CompVal).sum(axis=0))
            MCPerms[i,::,::]=(TS[shuffle_inds][:NumDraw]==CompVal).sum(axis=0)
        elif CompType=='lt':
            MCPerms[i]=(TS[shuffle_inds][:NumDraw]<CompVal).sum(axis=0)
        elif CompType=='le':
            MCPerms[i]=(TS[shuffle_inds][:NumDraw]<=CompVal).sum(axis=0)

    P=np.ones((I,J))*np.nan
    for i in range(I):
        for j in range(J):
            P[i,j]=stats.percentileofscore(MCPerms[::,i,j],obs[i,j])

    return P
