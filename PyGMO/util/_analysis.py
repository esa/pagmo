from PyGMO import *

class analysis:
    "This class will contain blahblah"

    def __init__(self,input_object):
        self.npoints=0
        self.points=[]
        self.f=[]
        self.grad_npoints=0
        self.grad_points=[]
        self.grad=[]
        self.c=[]
        self.local_nclusters=0
        self.local_initial_npoints=0

        object_type=str(type(input_object))
        if object_type=="<class 'PyGMO.core._core.population'>":
            self.prob=input_object.problem
            self.pop=input_object
        elif object_type[:21]=="<class 'PyGMO.problem":
            self.prob=input_object
            self.pop=[]
        else:
            raise ValueError(
             "analysis: input either a problem or a population object to initialise the class")


    def sample(self, npoints, method='sobol', first=1):
        self.points=[]
        self.f=[]
        self.npoints=npoints
        self.lb=list(self.prob.lb)
        self.ub=list(self.prob.ub)

        self.dim, self.cont_dim, self.int_dim, self.c_dim, self.ic_dim, self.f_dim = \
        self.prob.dimension, self.prob.dimension - self.prob.i_dimension, self.prob.i_dimension, self.prob.c_dimension, self.prob. ic_dimension, self.prob.f_dimension

        # if self.c_dim > 0:
        #     raise ValueError(
        #      "analysis.sample: this analyzer is not yet suitable for constrained optimisation")
        if self.npoints <= 0:
            raise ValueError(
             "analysis.sample: at least one point needs to be sampled")

        if method=='pop':
            poplength=len(self.pop)
            if poplength==0:
                raise ValueError(
                    "analysis.sample: method 'pop' specified but population object inexistant or void")
            elif poplength<npoints:
                raise ValueError(
                    "analysis.sample: it is not possible to sample more points than there are in the population via 'pop'")
            elif poplength==npoints:
                self.points=[list(self.pop[i].cur_x) for i in range(poplength)]
                self.f=[list(self.pop[i].cur_f) for i in range(poplength)]
            else:
                idx=range(poplength)
                try:
                    from numpy.random import randint
                except ImportError:
                    raise ImportError(
                        "analysis.sample needs numpy to run when sampling partially a population. Is it installed?")
                for i in range(poplength,poplength-npoints,-1):
                    r=idx.pop(randint(i))
                    self.points.append(list(self.pop[r].cur_x))
                    self.f.append(list(self.pop[r].cur_f))
        else:
            if method=='sobol':
                sampler=util.sobol(self.dim,first)
            elif method=='lhs':
                sampler=util.lhs(self.dim,npoints)
            else:
                raise ValueError(
                    "analysis.sample: method specified is not valid. choose 'sobol', 'lhs' or 'pop'")
            for i in range(npoints):
                temp=list(sampler.next()) #sample in the unit hypercube
                for j in range(self.dim):
                    temp[j]=temp[j]*self.ub[j]+(1-temp[j])*self.lb[j] #resize
                    if j>=self.cont_dim:
                        temp[j]=round(temp[j],0) #round if necessary
                self.points.append(temp)
                self.f.append(list(self.prob.objfun(temp)))

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

    def plot_f_distr(self):
        try:
            from scipy.stats import gaussian_kde
            from matplotlib.pyplot import plot,draw,title,show,legend
        except ImportError:
            raise ImportError(
                "analysis.plot_f_distr needs scipy and matplotlib to run. Are they installed?")
        if self.npoints==0:
            raise ValueError(
                "analysis.plot_f_distr: sampling first is necessary")
        plots=[]
        for i in range(self.f_dim):
            tmp=[]
            for j in range(self.npoints):
                tmp.append(self.f[j][i])
            x=sorted(tmp)
            kde=gaussian_kde(x)
            y=kde(x)
            plots.append(plot(x,y,label='objective '+str(i+1)))
        title('F-Distributions')
        legend()
        show(plots)

    def n_peaks_f(self,nf=0):
        try:
            from numpy import array,zeros
            from scipy.stats import gaussian_kde
        except ImportError:
            raise ImportError(
                "analysis.n_peaks_f needs scipy, numpy and matplotlib to run. Are they installed?")
        if self.npoints==0:
            raise ValueError(
                "analysis.n_peaks_f: sampling first is necessary")
        if nf==0:
            nf=self.npoints-1
        elif nf<3:
            raise ValueError(
                "analysis.n_peaks_f: value of nf too small")
        npeaks=[]
        for i in range(self.f_dim):
            npeaks.append(0)
            f=[a[i] for a in self.f]
            kde=gaussian_kde(f)
            df=self.ptp()[i]/nf
            x=[min(f)]
            for j in range(0,nf):
                x.append(x[j]+df)
            y=kde(x)
            minidx=[0]
            k=1
            for (a,b,c) in zip(y[0:nf-1],y[1:nf],y[2:nf+1]):
                if a>b<c:
                    minidx.append(k)
                k+=1
            minidx.append(nf+1)
            mode_mass=[kde.integrate_box_1d(x[minidx[j]],x[minidx[j+1]-1]) for j in range(len(minidx)-1)]
            for mode in mode_mass:
                if mode>0.01:
                    npeaks[i]+=1
        return npeaks

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
    def p_lin_conv(self,n_pairs=0,threshold=10**(-10),maximization=False):
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
            f_real=array(self.prob.objfun(x))
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
        self.lin_conv_npairs=n_pairs
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
            from numpy import corrcoef,transpose,dot
            from numpy.linalg import eigh
        except ImportError:
            raise ImportError(
                "analysis.f_correlation needs numpy to run. Is it installed?")
        M=corrcoef(self.f, rowvar=0)
        e=eigh(M)
        return (M.tolist(), e[0].tolist(), transpose(e[1]).tolist())
    
    def perform_f_pca(self,obj_corr=None,tc=0.95,tabs=0.1):
        try:
            from numpy import asarray,corrcoef,transpose,dot,argmax,argmin
            from numpy.linalg import eigh
            from itertools import combinations
        except ImportError:
            raise ImportError(
                "analysis.perform_f_pca needs numpy to run. Is it installed?")
        if obj_corr==None:
            obj_corr=self.f_correlation()
        M=obj_corr[0]
        eigenvals=obj_corr[1]
        eigenvects=obj_corr[2]
        #eigenvalue elimination of redundant objectives
        contributions=(asarray(eigenvals)/sum(eigenvals)).tolist()
        l=len(eigenvals)
        eig_order=[y for (x,y) in sorted(zip(contributions,range(l)),reverse=True)]
        cumulative_contribution=0
        keep=[]
        first=True
        for i in eig_order:
            index_p,index_n=argmax(eigenvects[i]),argmin(eigenvects[i])
            p,n=eigenvects[i][index_p],eigenvects[i][index_n]
            if first:
                first=False
                if p>0:
                    if all([k!=index_p for k in keep]):
                        keep.append(index_p)
                    if n<0:
                        if all([k!=index_n for k in keep]):
                            keep.append(index_n)
                else:
                    keep=range(l)
                    break
            elif eigenvals[i]<tabs:
                if abs(p)>abs(n):
                    if all([k!=index_p for k in keep]):
                        keep.append(index_p)
                else:
                    if all([k!=index_n for k in keep]):
                        keep.append(index_n)
            else:
                if n>=0:
                    if all([k!=index_p for k in keep]):
                        keep.append(index_p)
                elif p<=0:
                    keep=range(l)
                    break
                else:
                    if abs(n)>=p>=0.9*abs(n):
                        if all([k!=index_p for k in keep]):
                            keep.append(index_p)
                        if all([k!=index_n for k in keep]):    
                            keep.append(index_n)
                    elif p<0.9*abs(n):
                        if all([k!=index_n for k in keep]):
                            keep.append(index_n)
                    else:
                        if abs(n)>=0.8*p:
                            if all([k!=index_p for k in keep]):
                                keep.append(index_p)
                            if all([k!=index_n for k in keep]):    
                                keep.append(index_n)
                        else:
                            if all([k!=index_p for k in keep]):
                                keep.append(index_p)

            cumulative_contribution+=contributions[i]
            if cumulative_contribution>=tc or len(keep)==l:
                break
        #correlation elimination of redundant objectives
        if len(keep)>2: 
            c=list(combinations(keep,2))
            for i in range(len(c)):
                if all([x*y>0 for x,y in zip(M[c[i][0]],M[c[i][1]])]) and any([k==c[i][1] for k in keep]) and any([k==c[i][0] for k in keep]):
                    if keep.index(c[i][0])<keep.index(c[i][1]):
                        keep.remove(c[i][1])
                    else:
                        keep.remove(c[i][0])
        return sorted(keep)
      

#CURVATURE
#possible problem: tolerance needs to be relative to the magnitude of the result
    def get_gradient(self,sample_size=0,h=0.01,grad_tol=0.000001,zero_tol=0.000001):
        if self.npoints==0:
            raise ValueError(
                "analysis.get_gradient: sampling first is necessary")
        try:
            from numpy.random import randint
            from numpy import nanmean, asarray
        except ImportError:
            raise ImportError(
                "analysis.get_gradient needs numpy to run. Is it installed?")
        
        if sample_size<=0 or sample_size>=self.npoints:
            self.grad_points=range(self.npoints)
            self.grad_npoints=self.npoints
        else:
            self.grad_npoints=sample_size
            self.grad_points=[randint(self.npoints) for i in range(sample_size)] #avoid repetition?

        self.grad=[]
        self.grad_sparsity=0
        for i in self.grad_points:
            self.grad.append(self.richardson_gradient(x=self.points[i],h=h,grad_tol=grad_tol))
        self.average_abs_gradient=nanmean(abs(asarray(self.grad)),0)
        for i in range(self.f_dim):
            for j in range(self.dim):
                if abs(self.average_abs_gradient[i][j])<=zero_tol:
                    self.grad_sparsity+=1.
        self.grad_sparsity/=(self.dim*self.f_dim)

    def richardson_gradient(self,x,h,grad_tol,tmax=15):
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
                tmp=(array(self.prob.objfun(xu))-array(self.prob.objfun(xd)))/(2*hh)
                
                for j in range(self.f_dim):
                    d[t%2][0][j][i]=tmp[j]

            for k in range(1,t+1):
                d[t%2][k]=d[t%2][k-1]+(d[t%2][k-1]-d[(t+1)%2][k-1])/(4**k-1)

            if t>0:
                err=amax(abs(d[t%2][t]-d[(t+1)%2][t-1]))

            d[(t+1)%2].extend([zeros([self.f_dim,self.dim]),zeros([self.f_dim,self.dim])])
            t+=1

        return list(d[(t+1)%2][t-1])


    def get_hessian(self,sample_size=0,h=0.01,hess_tol=0.000001):
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
            self.hess.append(self.richardson_hessian(x=self.points[i],h=h,hess_tol=hess_tol))

    def richardson_hessian(self,x,h,hess_tol,tmax=15):
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

                    tmp=(array(self.prob.objfun(xu))-2*array(self.prob.objfun(x))+array(self.prob.objfun(xd)))/(hh**2)

                else:
                    xuu[ind[i][0]]+=hh
                    xuu[ind[i][1]]+=hh
                    xdd[ind[i][0]]-=hh
                    xdd[ind[i][1]]-=hh
                    xud[ind[i][0]]+=hh
                    xud[ind[i][1]]-=hh
                    xdu[ind[i][0]]-=hh
                    xdu[ind[i][1]]+=hh

                    tmp=(array(self.prob.objfun(xuu))-array(self.prob.objfun(xud))-array(self.prob.objfun(xdu))+array(self.prob.objfun(xdd)))/(4*hh*hh)

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

    def plot_gradient_sparsity(self,zero_tol=0.0001):
        if self.grad_npoints==0:
            raise ValueError(
                "analysis.plot_gradient_sparsity: sampling and getting gradient first is necessary")
        try:
            from matplotlib.pylab import spy,show,title,grid,xlabel,ylabel,xticks,yticks,draw
            from numpy import nanmean,asarray
        except ImportError:
            raise ImportError(
                "analysis.plot_gradient_sparsity needs matplotlib and numpy to run. Are they installed?")

        
        title('Gradient/Jacobian Sparsity ('+str(100*round(self.grad_sparsity,4))+'% sparse) \n \n')
        grid(True)
        xlabel('dimension')
        ylabel('objective')
        plot=spy(self.average_abs_gradient,precision=zero_tol,markersize=20)
        try:
            xlocs=range(self.dim)
            ylocs=range(self.f_dim)
            xlabels=[str(i) for i in range(1,self.dim+1)]
            ylabels=[str(i) for i in range(1,self.f_dim+1)]
            xticks(xlocs,[x.format(xlocs[i]) for i,x in enumerate(xlabels)])
            yticks(ylocs,[y.format(ylocs[i]) for i,y in enumerate(ylabels)])
        except IndexError, ValueError:
            pass
        show(plot)

    def plot_gradient_pcp(self,mode='x',scaled=True):
        if self.grad_npoints==0:
            raise ValueError(
                "analysis.plot_gradient_pcp: sampling and getting gradient first is necessary")
        if mode!='x' and mode!='f':
            raise ValueError(
                "analysis.plot_gradient_pcp: choose a valid value for mode ('x' or 'f')")
        if mode=='x' and self.dim==1:
            raise ValueError(
                "analysis.plot_gradient_pcp: mode 'x' makes no sense for univariate problems")
        if mode=='f' and self.f_dim==1:
            raise ValueError(
                "analysis.plot_gradient_pcp: mode 'f' makes no sense for single-objective problems")
        try:
            from pandas.tools.plotting import parallel_coordinates as pc
            from pandas import DataFrame as df
            from matplotlib.pyplot import show,title,grid,ylabel,xlabel
            from numpy import asarray,transpose
        except ImportError:
            raise ImportError(
                "analysis.plot_gradient_pcp needs pandas, numpy and matplotlib to run. Are they installed?")
        gradient=[]
        if scaled:
            ranges=self.ptp()
        if mode=='x':
            aux=0
            rowlabel=True
        else:
            aux=1
            rowlabel=False
        for i in range(self.grad_npoints):
            if rowlabel:
                tmp=[]
            else:
                tmp=[['x'+str(x+1) for x in range(self.dim)]]
            for j in range(self.f_dim):
                if rowlabel:
                    tmp.append(['objective '+str(j+1)])
                else:
                    tmp.append([])
                if scaled:
                    for k in range(self.dim):
                        tmp[j+aux].append(self.grad[i][j][k]*(self.ub[k]-self.lb[k])/ranges[j])
                else:
                    tmp[j+aux].extend(self.grad[i][j])
            if rowlabel:
                gradient.extend(tmp)
            else:
                tmp2=[]
                for ii in range(self.dim):
                    tmp2.append([])
                    for jj in range(self.f_dim+1):
                        tmp2[ii].append(tmp[jj][ii])
                gradient.extend(tmp2)
        gradient=df(gradient)

        title('Gradient PCP \n')
        grid(True)
        if scaled:
            scalelabel=' (scaled)'
        else:
            scalelabel=''
        ylabel('Derivative value'+scalelabel)
        if rowlabel:
            xlabel('Dimension')
        else:
            xlabel('Objective')
        plot=pc(gradient,0)
        show(plot)

#LOCAL SEARCH -> minimization assumed, single objective assumed

    def func(self,x,obj):

        return float(self.prob.objfun(x)[obj])

    def get_local_extrema0(self,sample_size=0,method='Powell'):
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
                res=minimize(self.func,self.points[i],(j,),method)
                tmp1.append(list(res.x))
                tmp2.append(res.nfev)
                tmp3.append(res.fun)
            self.local_extrema.append(tmp1)
            self.local_neval.append(tmp2)
            self.local_f.append(tmp3)

    def get_local_extrema(self,sample_size=0):#ADD METHOD CHOICE
        if self.npoints==0:
            raise ValueError(
                "analysis.get_local_extrema: sampling first is necessary")
        try:
            import numpy as np
        except ImportError:
            raise ImportError(
                "analysis.get_local_extrema needs numpy to run. Is it installed?")
        from numpy.random import randint
        from time import time
        
        if sample_size<=0 or sample_size>=self.npoints:
            self.local_initial_points=range(self.npoints)
            self.local_initial_npoints=self.npoints
        else:
            self.local_initial_npoints=sample_size
            self.local_initial_points=[randint(self.npoints) for i in range(sample_size)] #avoid repetition?

        self.local_extrema=[]
        #self.local_neval=[]// pygmo doesn't return it
        self.local_search_time=[]
        self.local_f=[]
        algo=algorithm.gsl_fr()
        if self.f_dim==1:
            decomposition=self.prob
        else:
            decomposition=problem.decompose(self.prob)
        for i in range(self.local_initial_npoints):
            pop=population(decomposition)
            pop.push_back(self.points[self.local_initial_points[i]])
            isl=island(algo,pop)
            start=time()
            isl.evolve(1)
            isl.join()
            finish=time()
            self.local_search_time.append((finish-start)*1000)
            self.local_extrema.append(isl.population.champion.x)
            self.local_f.append(isl.population.champion.f[0])

    def cluster_local_extrema(self,k=0,variance_ratio=0.99,single_cluster_tolerance=0.0001):
        if self.npoints==0:
            raise ValueError(
                "analysis_cluster_local_extrema: sampling first is necessary")
        if self.local_initial_npoints==0:
            raise ValueError(
                "analysis.cluster_local_extrema: getting local extrema first is necessary")
        try:
            import numpy as np
            import scipy as sp
            import sklearn as sk
        except ImportError:
            raise ImportError(
                "analysis.cluster_local_extrema needs numpy,scipy and sklearn to run. Are they installed?")
        from sklearn.cluster import KMeans

        dataset=np.zeros([self.local_initial_npoints,self.dim+1])#normalized dataset

        range_f=np.mean(np.ptp(self.f,0))
        mean_f=np.mean(self.f)

        if range_f<single_cluster_tolerance:
            raise ValueError(
                "analysis_cluster_local_extrema: the results appear to be constant")
        for i in range(self.local_initial_npoints):
            for j in range(self.dim):
                dataset[i][j]=(self.local_extrema[i][j]-0.5*self.ub[j]-0.5*self.lb[j])/(self.ub[j]-self.lb[j])
            dataset[i][self.dim]=(self.local_f[i]-mean_f)/range_f

        if k!=0:#cluster to given number of clusters
            clust=KMeans(k)

            #storage of output
            local_cluster=list(clust.fit_predict(dataset))
            self.local_nclusters=k
            cluster_size=np.zeros(k)
            for i in range(self.local_initial_npoints):
                cluster_size[local_cluster[i]]+=1
            cluster_size=list(cluster_size)

        else:#find out number of clusters
            clust=KMeans(1)
            total_distances=clust.fit_transform(dataset)
            total_center=clust.cluster_centers_[0]
            total_radius=max(total_distances)
            if total_radius<single_cluster_tolerance:#single cluster scenario
                #storage of output
                local_cluster=list(clust.predict(dataset))
                self.local_nclusters=1
                cluster_size=[0]
                for i in range(self.local_initial_npoints):
                    cluster_size[local_cluster[i]]+=1
                cluster_size=list(cluster_size)
            else:
                k=2 #multiple cluster scenario
                var_tot=sum([x**2 for x in total_distances])
                var_ratio=0
                while var_ratio<=variance_ratio:
                    clust=KMeans(k)
                    y=clust.fit_predict(dataset)
                    cluster_size=np.zeros(k)
                    var_exp=0
                    for i in range(self.local_initial_npoints):
                        cluster_size[y[i]]+=1
                    for i in range(k):
                        distance=np.linalg.norm(clust.cluster_centers_[i]-total_center)
                        var_exp+=cluster_size[i]*distance**2
                    var_ratio=var_exp/var_tot
                    k+=1
                #storage of output
                local_cluster=list(y)
                self.local_nclusters=k-1

        #more storage and reordering so clusters are ordered best to worst
        cluster_value=[clust.cluster_centers_[i][self.dim] for i in range(self.local_nclusters)]
        order=[x for (y,x) in sorted(zip(cluster_value,range(self.local_nclusters)))]

        self.local_cluster_x_centers=[]
        self.local_cluster_f_centers=[]
        self.local_cluster=[]
        self.local_cluster_size=[]
        for i in range(self.local_nclusters):
            self.local_cluster_size.append(cluster_size[order[i]])
            self.local_cluster_x_centers.append(clust.cluster_centers_[order[i]][:self.dim])
            self.local_cluster_f_centers.append(clust.cluster_centers_[order[i]][self.dim]*range_f+mean_f)
            for j in range(self.dim):
                self.local_cluster_x_centers[i][j]*=(self.ub[j]-self.lb[j])
                self.local_cluster_x_centers[i][j]+=0.5*(self.ub[j]+self.lb[j])
        for i in range(self.local_initial_npoints):
            for j in range(self.local_nclusters):
                if local_cluster[i]==order[j]:
                    self.local_cluster.append(j)
                    break

        #calculate cluster radius and center
        self.local_cluster_rx=[0 for i in range(self.local_nclusters)]
        f=[[] for i in range(self.local_nclusters)]
        for i in range(self.local_initial_npoints):
            c=self.local_cluster[i]
            if self.local_cluster_size[c]==1:
                f[c].append(0)
            else:
                rx=np.linalg.norm(np.asarray(self.local_extrema[i])-np.asarray(self.local_cluster_x_centers[c]))
                f[c].append(self.local_f[i])
                if rx>self.local_cluster_rx[c]:
                    self.local_cluster_rx[c]=rx
        self.local_cluster_df=[np.ptp(f[t],0).tolist() for t in range(self.local_nclusters)]

    def plot_local_cluster_pcp(self,together=True):
        if self.local_nclusters==0:
            raise ValueError(
                "analysis.plot_local_cluster_pcp: sampling, getting local extrema and clustering them first is necessary")
        if self.dim==1:
            raise ValueError(
                "analysis.plot_local_cluster_pcp: this makes no sense for univariate problems")
        try:
            from pandas.tools.plotting import parallel_coordinates as pc
            from pandas import DataFrame as df
            from matplotlib.pyplot import show,title,grid,ylabel,xlabel,legend,plot,subplot
            from numpy import asarray,transpose
        except ImportError:
            raise ImportError(
                "analysis.plot_gradient_pcp needs pandas, numpy and matplotlib to run. Are they installed?")
        if together:
            n=1
            dataset=[[[self.local_cluster[i]]+\
            [(self.points[self.local_initial_points[i]][j]-self.lb[j])/(self.ub[j]-self.lb[j]) for j in range(self.dim)]\
            for i in range(self.local_initial_npoints)]]
            dataset[0].sort()
            separatelabel=['' for i in range(self.local_nclusters)]
        else:
            n=self.local_nclusters
            dataset=[[] for i in range(self.local_nclusters)]
            for i in range(self.local_initial_npoints):
                dataset[self.local_cluster[i]].append([self.local_cluster[i]]+[(self.points[self.local_initial_points[i]][j]-self.lb[j])/(self.ub[j]-self.lb[j]) for j in range(self.dim)])
            separatelabel=[str(i) for i in range(self.local_nclusters)]
        for i in range(n):
            dataframe=df(dataset[i])
            title('Local Extrema Clusters PCP'+separatelabel[i]+' \n')
            grid(True)
            xlabel('Dimension')
            plot=pc(dataframe,0)
            show(plot)


#LEVELSET FEATURES (quite bad unless improved)
    def lda(self,threshold=50,tsp=0.1):
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
        from numpy.random import random
        clf=LDA()
        mce=[]
        for i in range (self.f_dim):
            per=self.percentile(threshold)[i]
            dataset=[[],[]]
            y=[[],[]]
            for j in range (self.npoints):
                r=random()
                if r<tsp:
                    index=1
                else:
                    index=0
                dataset[index].append(self.points[j])
                if self.f[j][i]>per:
                    y[index].append(1)
                else:
                    y[index].append(0)
            clf.fit(dataset[0],y[0])
            mce.append(1-clf.score(dataset[1],y[1]))
        return mce

    def qda(self,threshold=50,tsp=0.1):
        if self.npoints==0:
            raise ValueError(
                "analysis.qda: sampling first is necessary")
        try:
            import numpy as np
            import sklearn as sk
        except ImportError:
            raise ImportError(
                "analysis.qda needs numpy and scikit-learn to run. Are they installed?")
        from sklearn.qda import QDA
        from numpy import zeros
        from numpy.random import random
        clf=QDA()
        mce=[]
        for i in range (self.f_dim):
            per=self.percentile(threshold)[i]
            dataset=[[],[]]
            y=[[],[]]
            for j in range (self.npoints):
                r=random()
                if r<tsp:
                    index=1
                else:
                    index=0
                dataset[index].append(self.points[j])
                if self.f[j][i]>per:
                    y[index].append(1)
                else:
                    y[index].append(0)
            clf.fit(dataset[0],y[0])
            mce.append(1-clf.score(dataset[1],y[1]))
        return mce

    def kfdac(self,threshold=50,tsp=0.1):
        if self.npoints==0:
            raise ValueError(
                "analysis.kfdac: sampling first is necessary")
        try:
            import numpy as np
            import mlpy
        except ImportError:
            raise ImportError(
                "analysis.kfdac needs numpy and mlpy to run. Are they installed?")
        from mlpy import KFDAC, KernelGaussian
        from numpy import zeros
        from numpy.random import random
        K=KernelGaussian()
        clf=KFDAC(kernel=K)
        mce=[]
        for i in range (self.f_dim):
            per=self.percentile(threshold)[i]
            dataset=[[],[]]
            y=[[],[]]
            for j in range (self.npoints):
                r=random()
                if r<tsp:
                    index=1
                else:
                    index=0
                dataset[index].append(self.points[j])
                if self.f[j][i]>per:
                    y[index].append(1)
                else:
                    y[index].append(0)
            clf.learn(dataset[0],y[0])
            y_pred=clf.pred(dataset[1])
            score=(y_pred ==y[1])
            mce.append(1-np.mean(score))
        return mce

    def knn(self,threshold=50,tsp=0.1): #highly unuseful at the moment
        if self.npoints==0:
            raise ValueError(
                "analysis.knn: sampling first is necessary")
        try:
            import numpy as np
            import sklearn as sk
        except ImportError:
            raise ImportError(
                "analysis.knn needs numpy and scikit-learn to run. Are they installed?")
        from sklearn.neighbors import KNeighborsClassifier
        from numpy import zeros
        from numpy.random import random
        clf=KNeighborsClassifier(weights='distance')
        mce=[]
        for i in range (self.f_dim):
            per=self.percentile(threshold)[i]
            dataset=[[],[]]
            y=[[],[]]
            for j in range (self.npoints):
                r=random()
                if r<tsp:
                    index=1
                else:
                    index=0
                dataset[index].append(self.points[j])
                if self.f[j][i]>per:
                    y[index].append(1)
                else:
                    y[index].append(0)
            clf.fit(dataset[0],y[0])
            mce.append(1-clf.score(dataset[1],y[1]))
        return mce

#LEVELSET FEATURES WITH DAC, IMPROVED BUT IN PRINCIPLE WORSE THAN SVM
    def dac(self,threshold=50,classifier='k',k_test=10):
        if self.npoints==0:
            raise ValueError(
                "analysis.dac: sampling first is necessary")
        if classifier!='l' and classifier!= 'q' and classifier!='k':
            raise ValueError(
                "analysis.dac: choose a proper value for classifier ('l','q','k')")
        if threshold<=0 or threshold>=100:
            raise ValueError(
                "analysis.dac: threshold needs to be a value ]0,100[")
        try:
            import numpy as np
            import sklearn as sk
            import mlpy
        except ImportError:
            raise ImportError(
                "analysis.svm needs numpy, mlpy and scikit-learn to run. Are they installed?")
        from sklearn.cross_validation import StratifiedKFold, cross_val_score
        from sklearn.lda import LDA
        from sklearn.qda import QDA
        from sklearn.metrics import accuracy_score
        from mlpy import KFDAC, KernelGaussian

        if classifier=='l':
            clf=LDA()
        elif classifier=='q':
            clf=QDA()
        else:
            clf=KFDAC(kernel=KernelGaussian())
        per=self.percentile(threshold)
        
        dataset=[] #normalization of data
        for i in range(self.npoints):
            dataset.append(np.zeros(self.dim))
            for j in range(self.dim):
                dataset[i][j]=(self.points[i][j]-0.5*self.ub[j]-0.5*self.lb[j])/(self.ub[j]-self.lb[j])

        mce=[]
        for obj in range(self.f_dim):
            y=np.zeros(self.npoints) #classification of data
            for i in range(self.npoints):
                if self.f[i][obj]>per[obj]:
                    y[i]=1

            if classifier=='k':
                iterator=StratifiedKFold(y,k_test)
                i=0
                mce.append([])
                for train_index,test_index in iterator:
                    Xtrain=[]
                    Xtest=[]
                    ytrain=[]
                    ytest=[]
                    for i in train_index:
                        Xtrain.append(dataset[i])
                        ytrain.append(y[i])
                    for i in test_index:
                        Xtest.append(dataset[i])
                        ytest.append(y[i])
                    clf.learn(Xtrain,ytrain)
                    mce[obj].append(1-accuracy_score(clf.pred(Xtest),ytest))
            else:
                test_score=cross_val_score(estimator=clf,X=dataset,y=y,scoring=None,cv=StratifiedKFold(y,k_test))
                mce.append(list(np.ones(k_test)-test_score))

        return mce #mce[n_obj][k_test]

    def dac_p_values(self,threshold=50,k_test=10):
        if self.npoints==0:
            raise ValueError(
                "analysis.dac_p_values: sampling first is necessary")
        try:
            import scipy as sp
            import numpy as np
        except ImportError:
            raise ImportError(
                "analysis.dac_p_values needs scipy and numpy to run. Is it installed?")
        linear=self.dac(threshold=threshold, classifier='l', k_test=k_test)
        quadratic=self.dac(threshold=threshold, classifier='q', k_test=k_test)
        nonlinear=self.dac(threshold=threshold, classifier='k', k_test=k_test)
        l_q=[]
        q_n=[]
        l_n=[]
        for i in range(self.f_dim):
            l_q.append(sp.stats.mannwhitneyu(linear[i],quadratic[i])[1])
            l_n.append(sp.stats.mannwhitneyu(linear[i],nonlinear[i])[1])
            q_n.append(sp.stats.mannwhitneyu(quadratic[i],nonlinear[i])[1])
        return (list(np.mean(linear,1)),list(np.mean(quadratic,1)),list(np.mean(nonlinear,1)),l_q,l_n,q_n)

#LEVELSET FEATURES WITH SVM (PREFERABLY USE THESE...?)
    def svm(self,threshold=50,kernel='rbf',k_tune=3,k_test=10):
        if self.npoints==0:
            raise ValueError(
                "analysis.svm: sampling first is necessary")
        if kernel!='linear' and kernel!= 'quadratic' and kernel!='rbf':
            raise ValueError(
                "analysis.svm: choose a proper value for kernel ('linear','quadratic','rbf')")
        if threshold<=0 or threshold>=100:
            raise ValueError(
                "analysis.svm: threshold needs to be a value ]0,100[")
        try:
            import numpy as np
            import sklearn as sk
        except ImportError:
            raise ImportError(
                "analysis.svm needs numpy and scikit-learn to run. Are they installed?")
        from sklearn.cross_validation import StratifiedKFold, cross_val_score
        from sklearn.svm import SVC
        from sklearn.grid_search import GridSearchCV
        from sklearn.preprocessing import StandardScaler
        if kernel=='quadratic':
            kernel='poly'
        c_range=2.**np.arange(-5,16)
        if kernel=='linear':
            param_grid=dict(C=c_range)
        else:
            g_range=2.**np.arange(-15,4)
            param_grid=dict(gamma=g_range,C=c_range)
        per=self.percentile(threshold)
        
        dataset=[] #normalization of data
        for i in range(self.npoints):
            dataset.append(np.zeros(self.dim))
            for j in range(self.dim):
                dataset[i][j]=(self.points[i][j]-0.5*self.ub[j]-0.5*self.lb[j])/(self.ub[j]-self.lb[j])

        mce=[]
        for obj in range(self.f_dim):
            y=np.zeros(self.npoints) #classification of data
            for i in range(self.npoints):
                if self.f[i][obj]>per[obj]:
                    y[i]=1

            #grid search
            grid=GridSearchCV(estimator=SVC(kernel=kernel,degree=2),param_grid=param_grid,cv=StratifiedKFold(y,k_tune))
            grid.fit(dataset,y)
            #print grid.best_estimator_
            test_score=cross_val_score(estimator=grid.best_estimator_,X=dataset,y=y,scoring=None,cv=StratifiedKFold(y,k_test))
            mce.append(list(np.ones(k_test)-test_score))

        return mce #mce[n_obj][k_test]

    def svm_p_values(self,threshold=50,k_tune=3,k_test=10):
        if self.npoints==0:
            raise ValueError(
                "analysis.svm_p_values: sampling first is necessary")
        try:
            import scipy as sp
            import numpy as np
        except ImportError:
            raise ImportError(
                "analysis.svm_p_values needs scipy and numpy to run. Is it installed?")
        linear=self.svm(threshold=threshold, kernel='linear',k_tune=k_tune, k_test=k_test)
        quadratic=self.svm(threshold=threshold, kernel='quadratic',k_tune=k_tune, k_test=k_test)
        nonlinear=self.svm(threshold=threshold, kernel='rbf',k_tune=k_tune, k_test=k_test)
        l_q=[]
        q_n=[]
        l_n=[]
        for i in range(self.f_dim):
            l_q.append(sp.stats.mannwhitneyu(linear[i],quadratic[i])[1])
            l_n.append(sp.stats.mannwhitneyu(linear[i],nonlinear[i])[1])
            q_n.append(sp.stats.mannwhitneyu(quadratic[i],nonlinear[i])[1])
        return (list(np.mean(linear,1)),list(np.mean(quadratic,1)),list(np.mean(nonlinear,1)),l_q,l_n,q_n)


#CONSTRAINTS
    def compute_constraints(self):
        if self.npoints==0:
            raise ValueError(
                "analysis.compute_constraints: sampling first is necessary")
        self.c=[]
        if self.c_dim!=0:
            for i in range(self.npoints):
                self.c.append(list(self.prob.compute_constraints(self.points[i])))

    def ic_effectiveness(self):
        if self.npoints==0:
            raise ValueError(
                "analysis.constraint_feasibility: sampling first is necessary")
        ic_ef=[]
        if self.ic_dim!=0:
            if len(self.c)==0:
                raise ValueError(
                    "analysis.constraint_feasibility: compute constraints first")
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
        ec_f=[]
        if self.ic_dim-self.c_dim!=0:
            if len(self.c)==0:
                raise ValueError(
                    "analysis.constraint_feasibility: compute constraints first")
            for i in range(self.c_dim-self.ic_dim):
                ec_f.append(False)
                for j in range(self.npoints):
                    if self.c[j][i]==0 or (self.c[j][i]>0 and self.c[0][i]<0) or (self.c[0][i]>0 and self.c[j][i]<0):
                        ec_f[i]=True
            return ec_f

    def print_report(self,sample=100,p_f=[0,10,25,50,75,90,100],p_svm=[],p_dac=[],local_search=True): #UNFINISHED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        import numpy as np
        self.sample(sample)
        print "--------------------------------------------------------------------------------\n"
        print "PROBLEM PROPERTIES \n"
        print "--------------------------------------------------------------------------------\n"
        print "Dimension :                             ",self.dim," \n"
        print "     of which integer :                 ",self.int_dim," \n"
        print "Number of objectives :                  ",self.f_dim," \n"
        print "Number of constraints :                 ",self.c_dim," \n"
        print "     of which inequality constraints :  ",self.ic_dim," \n"
        print "Variable bounds : \n"
        for i in range(self.dim):
            print "     variable",i+1,":                       ","[",self.lb[i],",",self.ub[i],"]\n"
        print "Box constraints hypervolume :           ",self.box_hv()," \n"
        print "--------------------------------------------------------------------------------\n"
        print "F-DISTRIBUTION FEATURES (",self.f_dim," OBJECTIVES )\n"
        print "--------------------------------------------------------------------------------\n"
        print "Number of points sampled :              ",self.npoints," \n"
        print "Range of objective function :           ",list(self.ptp())," \n"
        print "Mean value :                            ",list(self.mean())," \n"
        print "Variance :                              ",list(self.var())," \n"
        if len(p_f)>0:
            print "Percentiles : \n"
            for i in p_f:
                if i<10:
                    print "     ",i,":                               ",list(self.percentile(i))," \n"
                elif i==100:
                    print "     ",i,":                             ",list(self.percentile(i))," \n"
                else:
                    print "     ",i,":                              ",list(self.percentile(i))," \n"  
        print "Skew :                                  ",list(self.skew())," \n"
        print "Kurtosis :                              ",list(self.kurtosis())," \n"
        print "Number of peaks of f-distribution :     ",self.n_peaks_f()," \n"
        print "--------------------------------------------------------------------------------\n"
        print "META-MODEL FEATURES \n"
        print "--------------------------------------------------------------------------------\n"
        print "Linear regression R2 :                  ",self.lin_reg()[1]," \n"
        if self.dim>=2:
            print "Linear regression with interaction R2 : ",self.lin_reg_inter()[1]," \n"
        print "Quadratic regression R2 :               ",self.poly_reg()[1]," \n"
        if len(p_svm)>0:
            print "--------------------------------------------------------------------------------\n"
            print "LEVELSET FEATURES : SVM\n"
            print "--------------------------------------------------------------------------------\n"
            for i in p_svm:
                svm_results=self.svm_p_values(threshold=i)
                print "Percentile",i,":\n"
                print "     Mean Misclassification Errors :\n"
                print "         Linear Kernel :                ",svm_results[0]," \n"
                print "         Quadratic Kernel :             ",svm_results[1]," \n"
                print "         Non-Linear Kernel (RBF):       ",svm_results[2]," \n"
                print "     P-Values : \n"
                print "         Linear/Quadratic :             ",svm_results[3]," \n"
                print "         Linear/Nonlinear :             ",svm_results[4]," \n"
                print "         Quadratic/Nonlinear :          ",svm_results[5]," \n"
        if len(p_dac)>0:
            print "--------------------------------------------------------------------------------\n"
            print "LEVELSET FEATURES : DAC\n"
            print "--------------------------------------------------------------------------------\n"
            for i in p_dac:
                dac_results=self.dac_p_values(threshold=i)
                print "Percentile",i,":\n"
                print "     Mean Misclassification Errors :\n"
                print "         LDA :                          ",dac_results[0]," \n"
                print "         QDA :                          ",dac_results[1]," \n"
                print "         KFDA (RBF Kernel):             ",dac_results[2]," \n"
                print "     P-Values : \n"
                print "         Linear/Quadratic :             ",dac_results[3]," \n"
                print "         Linear/Nonlinear :             ",dac_results[4]," \n"
                print "         Quadratic/Nonlinear :          ",dac_results[5]," \n"
        print "--------------------------------------------------------------------------------\n"
        print "PROBABILITY OF LINEARITY AND CONVEXITY\n"
        print "--------------------------------------------------------------------------------\n"
        p=self.p_lin_conv()
        print "Number of pairs of points used :        ",self.lin_conv_npairs," \n"
        print "Probability of linearity :              ",p[0]," \n"
        print "Probability of convexity :              ",p[1]," \n"
        print "Mean deviation from linearity :         ",p[2]," \n"
        if self.c_dim==0 and local_search:
            print "--------------------------------------------------------------------------------\n"
            print "LOCAL SEARCH\n"
            print "--------------------------------------------------------------------------------\n"
            self.get_local_extrema()
            self.cluster_local_extrema()
            print "Local searches performed :              ",self.local_initial_npoints," \n"
            print "Quartiles of CPU time per search [ms]:  ",round(np.percentile(self.local_search_time,0),3),"/",round(np.percentile(self.local_search_time,25),3),"/",round(np.percentile(self.local_search_time,50),3),"/",round(np.percentile(self.local_search_time,75),3),"/",round(np.percentile(self.local_search_time,100),3)," \n"
            print "Number of clusters identified :         ",self.local_nclusters," \n"
            print "Cluster properties (max. best 5 clusters): \n"
            for i in range(min((self.local_nclusters,5))):
                print "     Cluster n.:                        ",i+1," \n"
                print "         Size:                          ",self.local_cluster_size[i],", ",100*round(self.local_cluster_size[i]/self.local_initial_npoints,4),"% \n"
                print "         Mean objective value :         ",self.local_cluster_f_centers[i]," \n"
                print "         Cluster center :               ",self.local_cluster_x_centers[i]," \n"
                print "         Cluster diameter in F :        ",self.local_cluster_df[i]," \n"
                print "         Cluster radius in X :          ",self.local_cluster_rx[i]," \n"
        print "--------------------------------------------------------------------------------\n"
        print "CURVATURE : GRADIENT/JACOBIAN \n"
        print "--------------------------------------------------------------------------------\n"
        self.get_gradient()
        print "Number of points evaluated :            ",self.grad_npoints," \n"
        print "Gradient sparsity :                     ",100*round(self.grad_sparsity,4),"% \n"
        print "--------------------------------------------------------------------------------\n"
        print "CURVATURE : HESSIAN \n"
        print "--------------------------------------------------------------------------------\n"
        self.get_hessian()
        print "Number of points evaluated :            ",self.hess_npoints," \n"
        
        if self.f_dim>1:
            print "--------------------------------------------------------------------------------\n"
            print "OBJECTIVE CORRELATION \n"
            print "--------------------------------------------------------------------------------\n"
            obj_corr=self.f_correlation()
            critical_obj=self.perform_f_pca(obj_corr)
            print "Objective correlation matrix :          ",obj_corr[0][0]," \n"
            for i in range(1,self.f_dim):
                print "                                        ",obj_corr[0][i]," \n"
            print "Eigenvalues (squared) :                 ",obj_corr[1]," \n"
            print "Eigenvalue relative contribution :      ",[str(100*round(i,4))+'%' for i in (np.asarray(obj_corr[1])/sum(obj_corr[1])).tolist()]," \n"
            print "Eigenvectors :                          ",obj_corr[2][0]," \n"
            for i in range(1,self.f_dim):
                print "                                        ",obj_corr[2][i]," \n"
            print "Critical objectives from first PCA :    ",[int(i) for i in (np.ones(len(critical_obj))+np.asarray(critical_obj)).tolist()]," \n"
        if self.c_dim>0:
            print "--------------------------------------------------------------------------------\n"
            print "CONSTRAINT EFFECTIVENESS/FEASIBILITY \n"
            print "--------------------------------------------------------------------------------\n"
            self.compute_constraints()
            ic_ef=self.ic_effectiveness()
            ec_f=self.ec_feasibility()
            if self.c_dim!=self.ic_dim:
                print "     Equality constraint      |      Feasibility   \n"
            for i in range(self.c_dim-self.ic_dim):
                print "            ",i+1,"                        ",ec_f[i],"       \n"
            if self.ic_dim>0:
                print "    Inequality constraint     |     Effectiveness   \n"
            for i in range(self.ic_dim):
                print "            ",self.c_dim-self.ic_dim+i+1,"                       ",100*round(ic_ef[i],4),"%      \n"
        #PLOTS 
        self.plot_f_distr()
        self.plot_gradient_sparsity()
        if self.dim>1:
            self.plot_gradient_pcp('x')
        if self.f_dim>1:
            self.plot_gradient_pcp('f')