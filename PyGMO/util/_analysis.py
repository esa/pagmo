
from PyGMO import *

class analysis:
    "This class will contain blahblah"

    def __init__(self):
            self.npoints=0
            self.points=[]
            self.f=[]
            self.grad_npoints=0
            self.grad_points=[]
            self.grad=[]
            self.c=[]


    def sample(self, prob, npoints, method='sobol', first=1):
        self.points=[]
        self.f=[]
        self.npoints=npoints
        self.lb=list(prob.lb)
        self.ub=list(prob.ub)

        self.dim, self.cont_dim, self.int_dim, self.c_dim, self.ic_dim, self.f_dim = \
        prob.dimension, prob.dimension - prob.i_dimension, prob.i_dimension, prob.c_dimension, prob. ic_dimension, prob.f_dimension

        # if self.c_dim > 0:
        #     raise ValueError(
        #      "analysis.sample: this analyzer is not yet suitable for constrained optimisation")
        if self.npoints == 0:
            raise ValueError(
             "analysis.sample: at least one point needs to be sampled")

        if method=='sobol':
            sampler=util.sobol(self.dim,first)
        elif method=='lhs':
            sampler=util.lhs(self.dim,npoints)
        else:
            raise ValueError(
                "analysis.sample: method specified is not valid. choose 'sobol' or 'lhs'")

        for i in range(npoints):
            temp=list(sampler.next()) #sample in the unit hypercube
            for j in range(self.dim):
                temp[j]=temp[j]*self.ub[j]+(1-temp[j])*self.lb[j] #resize
                if j>=self.cont_dim:
                    temp[j]=round(temp[j],0) #round if necessary
            self.points.append(temp)
            self.f.append(list(prob.objfun(temp)))

#f-DISTRIBUTION FEATURES
   
    def skew(self):
        try:
            import scipy as sp
        except ImportError:
            raise ImportError(
                "analysis.skew needs scipy to run. Is it installed?")
        if self.npoints==0:
            raise ValueError(
                "analysis.skew: sampling first is necessary")
        from scipy.stats import skew
        return skew(self.f)

    def kurtosis(self):
        try:
            import scipy as sp
        except ImportError:
            raise ImportError(
                "analysis.kurtosis needs scipy to run. Is it installed?")
        if self.npoints==0:
            raise ValueError(
                "analysis.kurtosis: sampling first is necessary")
        from scipy.stats import kurtosis
        return kurtosis(self.f)

    def mean(self):
        try:
            import numpy as np
        except ImportError:
            raise ImportError(
                "analysis.mean needs numpy to run. Is it installed?")
        if self.npoints==0:
            raise ValueError(
                "analysis.mean: sampling first is necessary")
        from numpy import mean
        return mean(self.f,0)

    def var(self):
        try:
            import numpy as np
        except ImportError:
            raise ImportError(
                "analysis.var needs numpy to run. Is it installed?")
        if self.npoints==0:
            raise ValueError(
                "analysis.var: sampling first is necessary")
        from numpy import var
        return var(self.f,0)

    def ptp(self):
        try:
            import numpy as np
        except ImportError:
            raise ImportError(
                "analysis.ptp needs numpy to run. Is it installed?")
        if self.npoints==0:
            raise ValueError(
                "analysis.ptp: sampling first is necessary")
        from numpy import ptp
        return ptp(self.f,0)

    def percentile(self,p):
        try:
            import numpy as np
        except ImportError:
            raise ImportError(
                "analysis.percentile needs numpy to run. Is it installed?")
        if self.npoints==0:
            raise ValueError(
                "analysis.percentile: sampling first is necessary")
        from numpy import percentile
        return percentile(self.f,p,0)

    #ADD MORE FUNCTIONS IF NECESSARY

#BOX CONSTRAINT HYPERVOLUME COMPUTATION
    def box_hv(self):
        if self.npoints==0:
            raise ValueError(
                "analysis.box_hv: sampling first is necessary")
        hv=1
        for i in range(self.dim):
            hv=hv*(self.ub[i]-self.lb[i])
        return hv

#DEGREE OF LINEARITY AND CONVEXITY
    def p_lin_conv(self,prob,n_pairs=0,threshold=10**(-10),maximization=False):
        if self.npoints==0:
            raise ValueError(
                "analysis.p_lin_conv: sampling first is necessary")
        if n_pairs==0:
            n_pairs=self.npoints
        try:
            import numpy as np
        except ImportError:
            raise ImportError(
                "analysis.p_lin_conv needs numpy to run. Is it installed?")

        from numpy.random import random,randint
        from numpy import array
        p_lin=np.zeros(self.f_dim)
        p_conv=np.zeros(self.f_dim)
        mean_dev=np.zeros(self.f_dim)
        for i in range(n_pairs):
            i1=randint(self.npoints)
            i2=randint(self.npoints)
            while (i2==i1):
                i2=randint(self.npoints)
            r=random()
            x=r*array(self.points[i1])+(1-r)*array(self.points[i2])
            f_lin=r*array(self.f[i1])+(1-r)*array(self.f[i2])
            f_real=array(prob.objfun(x))
            delta=f_lin-f_real
            mean_dev+=abs(delta)
            if maximization:
                delta*=-1
            for j in range(self.f_dim):
                if abs(delta[j])<threshold:
                    p_lin[j]+=1
                elif delta[j]>0:
                    p_conv[j]+=1
        p_lin/=n_pairs
        p_conv/=n_pairs
        mean_dev/=n_pairs
        return (list(p_lin),list(p_conv),list(mean_dev))

#META-MODEL FEATURES
    def lin_reg(self):
        if self.npoints==0:
            raise ValueError(
                "analysis.lin_reg: sampling first is necessary")
        try:
            import numpy as np
        except ImportError:
            raise ImportError(
                "analysis.lin_reg needs numpy to run. Is it installed?")
        from numpy import array
        from numpy.linalg import lstsq
        A=[]
        w=[]
        sst=self.var()*self.npoints
        m=self.mean()
        ssr=np.zeros(self.f_dim)
        for i in range(self.npoints):
            A.append(self.points[i]+[1])
        A=array(A)
        for i in range(self.f_dim):
            b=[]
            for j in range(self.npoints):
                b.append(self.f[j][i])
            b=array(b)
            temp=lstsq(A,b)[0]
            w.append(list(temp))
            for j in range(self.npoints):
                ssr[i]+=(np.dot(temp,A[j])-m[i])**2
        r2=list(ssr/sst)
        return (w,r2)


    def lin_reg_inter(self,interaction_order=2):
        if self.npoints==0:
            raise ValueError(
                "analysis.lin_reg_corr: sampling first is necessary")
        if interaction_order<2 or interaction_order>self.dim:
            raise ValueError(
                "analysis.lin_reg_corr: interaction order should be in range [2,dim]")
        try:
            import numpy as np
        except ImportError:
            raise ImportError(
                "analysis.lin_reg_corr needs numpy to run. Is it installed?")
        from numpy import array
        from numpy.linalg import lstsq
        from itertools import combinations
        A=[]
        w=[]
        inter=[]
        for i in range(interaction_order,1,-1):
            inter=inter+list(combinations(range(self.dim),i))
        n_inter=len(inter)
        sst=self.var()*self.npoints
        m=self.mean()
        ssr=np.zeros(self.f_dim)
        for i in range(self.npoints):
            c=[]
            for j in range(n_inter):
                prod=1
                for k in range(len(inter[j])):
                    prod*=self.points[i][inter[j][k]]
                c.append(prod)
            A.append(c+self.points[i]+[1])
        A=array(A)
        for i in range(self.f_dim):
            b=[]
            for j in range(self.npoints):
                b.append(self.f[j][i])
            b=array(b)
            temp=lstsq(A,b)[0]
            w.append(list(temp))
            for j in range(self.npoints):
                ssr[i]+=(np.dot(temp,A[j])-m[i])**2
        r2=list(ssr/sst)
        return (w,r2)

    def poly_reg(self,regression_degree=2):
        if self.npoints==0:
            raise ValueError(
                "analysis.lin_reg_corr: sampling first is necessary")
        try:
            import numpy as np
        except ImportError:
            raise ImportError(
                "analysis.lin_reg_corr needs numpy to run. Is it installed?")
        from numpy import array
        from numpy.linalg import lstsq
        from itertools import combinations_with_replacement
        A=[]
        w=[]
        coef=[]
        for i in range(regression_degree,1,-1):
            coef=coef+list(combinations_with_replacement(range(self.dim),i))
        n_coef=len(coef)
        sst=self.var()*self.npoints
        m=self.mean()
        ssr=np.zeros(self.f_dim)
        for i in range(self.npoints):
            c=[]
            for j in range(n_coef):
                prod=1
                for k in range(len(coef[j])):
                    prod*=self.points[i][coef[j][k]]
                c.append(prod)
            A.append(c+self.points[i]+[1])
        A=array(A)
        for i in range(self.f_dim):
            b=[]
            for j in range(self.npoints):
                b.append(self.f[j][i])
            b=array(b)
            temp=lstsq(A,b)[0]
            w.append(list(temp))
            for j in range(self.npoints):
                ssr[i]+=(np.dot(temp,A[j])-m[i])**2
        r2=list(ssr/sst)
        return (w,r2)

#OBJECTIVE REDUNDANCY
    def f_correlation(self):
        if self.npoints==0:
            raise ValueError(
                "analysis.f_correlation: sampling first is necessary")
        if self.f_dim<2:
            raise ValueError(
                "analysis.f_correlation: this test makes no sense for single-objective optimisation")
        try:
            import numpy as np
        except ImportError:
            raise ImportError(
                "analysis.f_correlation needs numpy to run. Is it installed?")
        from numpy import corrcoef
        M=corrcoef(self.f, rowvar=0)
        return (M, np.linalg.eigh(M)[0], list(np.transpose(np.linalg.eigh(M)[1])))
    
    #enhance to perform pca...!?

#CURVATURE
#problem: tolerance needs to be relative to the magnitude of the result
    def get_gradient(self,prob,sample_size=0,h=0.01,grad_tol=0.000001,zero_tol=0.0001):
        if self.npoints==0:
            raise ValueError(
                "analysis.get_gradient: sampling first is necessary")
        try:
            import numpy as np
        except ImportError:
            raise ImportError(
                "analysis.get_gradient needs numpy to run. Is it installed?")
        from numpy.random import randint
        
        if sample_size<=0 or sample_size>=self.npoints:
            self.grad_points=range(self.npoints)
            self.grad_npoints=self.npoints
        else:
            self.grad_npoints=sample_size
            self.grad_points=[randint(self.npoints) for i in range(sample_size)] #avoid repetition?

        self.grad=[]
        self.grad_sparsity=np.zeros(self.f_dim)
        for i in self.grad_points:
            self.grad.append(self.richardson_gradient(x=self.points[i],prob=prob,h=h,grad_tol=grad_tol))
            for j in range(self.f_dim):
                for k in range(self.dim):
                    if abs(self.grad[i][j][k])<=zero_tol:
                        self.grad_sparsity[j]+=1.
        self.grad_sparsity/=(self.grad_npoints*self.dim*self.f_dim)

    def richardson_gradient(self,x,prob,h,grad_tol,tmax=15):
        from numpy import array, zeros, amax
        d=[[zeros([self.f_dim,self.dim])],[]]
        hh=2*h
        err=1
        t=0
        while (err>grad_tol and t<tmax):
            hh/=2
            for i in range(self.dim):
                xu=list(x)
                xd=list(x)
                xu[i]+=hh
                xd[i]-=hh
                tmp=(array(prob.objfun(xu))-array(prob.objfun(xd)))/(2*hh)
                
                for j in range(self.f_dim):
                    d[t%2][0][j][i]=tmp[j]

            for k in range(1,t+1):
                d[t%2][k]=d[t%2][k-1]+(d[t%2][k-1]-d[(t+1)%2][k-1])/(4**k-1)

            if t>0:
                err=amax(abs(d[t%2][t]-d[(t+1)%2][t-1]))

            d[(t+1)%2].extend([zeros([self.f_dim,self.dim]),zeros([self.f_dim,self.dim])])
            t+=1

        return list(d[(t+1)%2][t-1])


    def get_hessian(self,prob,sample_size=0,h=0.01,hess_tol=0.000001,zero_tol=0.0001):
        if self.npoints==0:
            raise ValueError(
                "analysis.get_hessian: sampling first is necessary")
        try:
            import numpy as np
        except ImportError:
            raise ImportError(
                "analysis.get_hessian needs numpy to run. Is it installed?")
        from numpy.random import randint
        
        if sample_size<=0 or sample_size>=self.npoints:
            self.hess_points=range(self.npoints)
            self.hess_npoints=self.npoints

        else:
            self.hess_npoints=sample_size
            self.hess_points=[randint(self.npoints) for i in range(sample_size)] #avoid repetition?

        self.hess=[]
        for i in self.hess_points:
            self.hess.append(self.richardson_hessian(x=self.points[i],prob=prob,h=h,hess_tol=hess_tol))

    def richardson_hessian(self,x,prob,h,hess_tol,tmax=15):
        from numpy import array, zeros, amax
        from itertools import combinations_with_replacement
        ind=list(combinations_with_replacement(range(self.dim),2))
        n_ind=len(ind)
        d=[[zeros([self.f_dim,n_ind])],[]]
        hh=2*h
        err=1
        t=0
        while (err>hess_tol and t<tmax):
            hh/=2
            for i in range(n_ind):
                xu=list(x)
                xd=list(x)
                xuu=list(x)
                xdd=list(x)
                xud=list(x)
                xdu=list(x)

                if ind[i][0]==ind[i][1]:
                    xu[ind[i][0]]+=hh
                    xd[ind[i][0]]-=hh

                    tmp=(array(prob.objfun(xu))-2*array(prob.objfun(x))+array(prob.objfun(xd)))/(hh**2)

                else:
                    xuu[ind[i][0]]+=hh
                    xuu[ind[i][1]]+=hh
                    xdd[ind[i][0]]-=hh
                    xdd[ind[i][1]]-=hh
                    xud[ind[i][0]]+=hh
                    xud[ind[i][1]]-=hh
                    xdu[ind[i][0]]-=hh
                    xdu[ind[i][1]]+=hh

                    tmp=(array(prob.objfun(xuu))-array(prob.objfun(xud))-array(prob.objfun(xdu))+array(prob.objfun(xdd)))/(4*hh*hh)

                for j in range(self.f_dim):
                    d[t%2][0][j][i]=tmp[j]

            for k in range(1,t+1):
                d[t%2][k]=d[t%2][k-1]+(d[t%2][k-1]-d[(t+1)%2][k-1])/(4**k-1)

            if t>0:
                err=amax(abs(d[t%2][t]-d[(t+1)%2][t-1]))

            d[(t+1)%2].extend([zeros([self.f_dim,n_ind]),zeros([self.f_dim,n_ind])])
            t+=1

        hessian=[]
        for i in range(self.f_dim):
            mat=zeros([self.dim,self.dim])
            for j in range(n_ind):
                mat[ind[j][0]][ind[j][1]]=d[(t+1)%2][t-1][i][j]
                mat[ind[j][1]][ind[j][0]]=d[(t+1)%2][t-1][i][j]
            hessian.append(mat)

        return hessian

#LOCAL SEARCH -> minimization assumed, single objective assumed
    def func(self,x,prob,obj):
        return float(prob.objfun(x)[obj])

    def get_local_extrema(self,prob,sample_size=0,method='Powell'):
        if self.npoints==0:
            raise ValueError(
                "analysis.get_local_extrema: sampling first is necessary")
        if (method=='Powell' or method=='Nelder-Mead' or method=='BFGS' or method=='CG')==False:
            raise ValueError(
                "analysis.get_local_extrema: choose a method amongst 'Powell', 'Nelder-Mead', 'BFGS' or 'CG'")
        try:
            import scipy as sp
            import numpy as np
        except ImportError:
            raise ImportError(
                "analysis.get_local_extrema needs numpy and scipy to run. Are they installed?")
        from numpy.random import randint
        from scipy.optimize import minimize
        
        if sample_size<=0 or sample_size>=self.npoints:
            self.local_initial_points=range(self.npoints)
            self.local_initial_npoints=self.npoints
        else:
            self.local_initial_npoints=sample_size
            self.local_initial_points=[randint(self.npoints) for i in range(sample_size)] #avoid repetition?

        self.local_extrema=[]
        self.local_neval=[]
        self.local_f=[]
        for j in range(self.f_dim):
            tmp1=[]
            tmp2=[]
            tmp3=[]
            for i in self.local_initial_points:
                res=minimize(self.func,self.points[i],(prob,j),method)
                tmp1.append(list(res.x))
                tmp2.append(res.nfev)
                tmp3.append(res.fun)
            self.local_extrema.append(tmp1)
            self.local_neval.append(tmp2)
            self.local_f.append(tmp3)

#LEVELSET FEATURES
    def lda(self,threshold=50):
        if self.npoints==0:
            raise ValueError(
                "analysis.lda: sampling first is necessary")
        try:
            import numpy as np
            import sklearn as sk
        except ImportError:
            raise ImportError(
                "analysis.lda needs numpy and scikit-learn to run. Are they installed?")
        from sklearn.lda import LDA
        from numpy import zeros
        clf=LDA()
        mce=[]
        for i in range (self.f_dim):
            y=zeros(self.npoints)
            for j in range (self.npoints):
                if self.f[j][i]>self.percentile(threshold)[i]:
                    y[j]=1
            clf.fit(self.points,y)
            mce.append(1-clf.score(self.points,y))
        return mce

    def qda(self,threshold=50):
        if self.npoints==0:
            raise ValueError(
                "analysis.lda: sampling first is necessary")
        try:
            import numpy as np
            import sklearn as sk
        except ImportError:
            raise ImportError(
                "analysis.lda needs numpy and scikit-learn to run. Are they installed?")
        from sklearn.qda import QDA
        from numpy import zeros
        clf=QDA()
        mce=[]
        for i in range (self.f_dim):
            y=zeros(self.npoints)
            for j in range (self.npoints):
                if self.f[j][i]>self.percentile(threshold)[i]:
                    y[j]=1
            clf.fit(self.points,y)
            mce.append(1-clf.score(self.points,y))
        return mce

#CONSTRAINTS
    def compute_constraints(self,prob):
        if self.npoints==0:
            raise ValueError(
                "analysis.compute_constraints: sampling first is necessary")
        self.c=[]
        if self.c_dim!=0:
            for i in range(self.npoints):
                self.c.append(list(prob.compute_constraints(self.points[i])))

    def ic_effectiveness(self):
        if self.npoints==0:
            raise ValueError(
                "analysis.constraint_feasibility: sampling first is necessary")
        if self.ic_dim!=0:
            if len(self.c)==0:
                raise ValueError(
                    "analysis.constraint_feasibility: compute constraints first")
            ic_ef=[]
            for i in range(self.ic_dim):
                ic_ef.append(0)
            dp=1./self.npoints
            for i in range(self.npoints):
                for j in range(-self.ic_dim,0):
                    if self.c[i][j]>=0:
                        ic_ef[j]+=dp
            return ic_ef

    def ec_feasibility(self):
        if self.npoints==0:
            raise ValueError(
                "analysis.constraint_feasibility: sampling first is necessary")
        if self.ic_dim-self.c_dim!=0:
            if len(self.c)==0:
                raise ValueError(
                    "analysis.constraint_feasibility: compute constraints first")
            ec_f=[]
            for i in range(self.c_dim-self.ic_dim):
                ec_f.append(1)
                for j in range(self.npoints):
                    if self.c[j][i]==0 or (self.c[j][i]>0 and self.c[0][i]<0) or (self.c[0][i]>0 and self.c[j][i]<0):
                        ec_f[i]=0
            return ec_f


    def print_report(self,prob): #UNFINISHED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        print "---------------------------------------------------------------------\n"
        print "PROBLEM PROPERTIES \n"
        print "---------------------------------------------------------------------\n"
        print "Dimension :                             ",self.dim," \n"
        print "     continuous :                       ",self.cont_dim," \n"
        print "     integer :                          ",self.int_dim," \n"
        print "Number of objectives :                  ",self.f_dim," \n"
        print "Number of constraints :                 ",self.c_dim," \n"
        print "     of which inequality constraints :  ",self.ic_dim," \n"
        print "Bounds : \n"
        for i in range(self.dim):
            print "     variable",i,":                         ","[",self.lb[i],",",self.ub[i],"]\n"
        #print "Box constraints hypervolume :           ",self.box_hv()," \n"
        print "---------------------------------------------------------------------\n"
        print "F-DISTRIBUTION FEATURES (",self.f_dim," OBJECTIVES )\n"
        print "---------------------------------------------------------------------\n"
        print "Number of points sampled :              ",self.npoints," \n \n"
        print "Range :                                 ",list(self.ptp())," \n"
        print "Mean value :                            ",list(self.mean())," \n"
        print "Variance :                              ",list(self.var())," \n"
        print "Percentiles :"
        print "     0 :                                ",list(self.percentile(0))," \n"
        print "    10 :                                ",list(self.percentile(10))," \n"
        print "    25 :                                ",list(self.percentile(25))," \n"
        print "    50 :                                ",list(self.percentile(50))," \n"
        print "    75 :                                ",list(self.percentile(75))," \n"
        print "    90 :                                ",list(self.percentile(90))," \n"
        print "   100 :                                ",list(self.percentile(100))," \n"
        print "Skew :                                  ",list(self.skew())," \n"
        print "Kurtosis :                              ",list(self.kurtosis())," \n"
        print "---------------------------------------------------------------------\n"
        print "META-MODEL FEATURES \n"
        print "---------------------------------------------------------------------\n"
        print "Linear regression R2 :                  ",self.lin_reg()[1]," \n"
        print "Linear regression with interaction R2 : ",self.lin_reg_inter()[1]," \n"
        print "Quadratic regression R2 :               ",self.poly_reg()[1]," \n"
        print "---------------------------------------------------------------------\n"
        print "LEVELSET FEATURES \n"
        print "---------------------------------------------------------------------\n"
        print "LDA Misclassification error : \n"
        print "     Percentile 10 :                    ",self.lda(10)," \n"
        print "     Percentile 25 :                    ",self.lda(25)," \n"
        print "     Percentile 50 :                    ",self.lda(50)," \n"
        print "QDA Misclassification error : \n"
        print "     Percentile 10 :                    ",self.qda(10)," \n"
        print "     Percentile 25 :                    ",self.qda(25)," \n"
        print "     Percentile 50 :                    ",self.qda(50)," \n"
        print "---------------------------------------------------------------------\n"
        print "PROBABILITY OF LINEARITY AND CONVEXITY\n"
        print "---------------------------------------------------------------------\n"
        p=self.p_lin_conv(prob)
        print "Probability of linearity :              ",p[0]," \n"
        print "Probability of convexity :              ",p[1]," \n"
        print "Mean deviation from linearity :         ",p[2]," \n"
