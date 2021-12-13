"""A Python class containing an implimentation of MPPCA denoising.
        
    Inputs are a 4D image with dimentions (X x Y x Z x N)
    
    Usage:
    import mpdenoise as mp
    imgdn, sigma, nparameters = mp.denoise(img, kernel=[7,7,7], M=60, shrinkage='frobenius')
    
    LICENCE
    Authors: Benjamin Ades-Aron (Benjamin.Ades-Aron@nyulangone.org)
    Copyright (c) 2016 New York University
    
    Permission is hereby granted, free of charge, to any non-commercial entity
    ('Recipient') obtaining a copy of this software and associated
    documentation files (the 'Software'), to the Software solely for
    non-commercial research, including the rights to use, copy and modify the
    Software, subject to the following conditions:
    
    1. The above copyright notice and this permission notice shall be
    included by Recipient in all copies or substantial portions of the
    Software.
    
    2. THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIESOF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
    NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BELIABLE FOR ANY CLAIM,
    DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
    OTHERWISE, ARISING FROM, OUT OF ORIN CONNECTION WITH THE SOFTWARE OR THE
    USE OR OTHER DEALINGS IN THE SOFTWARE.
    
    3. In no event shall NYU be liable for direct, indirect, special,
    incidental or consequential damages in connection with the Software.
    Recipient will defend, indemnify and hold NYU harmless from any claims or
    liability resulting from the use of the Software by recipient.

    4. Neither anything contained herein nor the delivery of the Software to
    recipient shall be deemed to grant the Recipient any right or licenses
    under any patents or patent application owned by NYU.

    5. The Software may only be used for non-commercial research and may not
    be used for clinical care.

    6. Any publication by Recipient of research involving the Software shall
    cite the references listed below.

    REFERENCES
    Veraart, J.; Fieremans, E. & Novikov, D.S. Diffusion MRI noise mapping
    using random matrix theory Magn. Res. Med., 2016, early view, doi:
    10.1002/mrm.26059
    """

import os, sys
import numpy as np
import multiprocessing
from joblib import Parallel, delayed
from tqdm import tqdm
import scipy.linalg
#import jax.scipy as jsp
#from jax import jit

class MP(object):
    def __init__(self, dwi, kernel, patchtype, patchsize, shrinkage, algorithm, rotation, pmap, crop):
        if dwi.ndim > 4:
            self.coil = True
        else:
            self.coil = False

        if not rotation:
            self.rotation = False
        else:
            self.rotation = True

        if shrinkage is None:
            self.shrink = 'threshold'
        else:
            self.shrink = shrinkage

        if algorithm is None:
            self.algo = 'jespersen'
        else:
            self.algo = algorithm

        if crop is None:
            self.tn = 5
            crop = 0
        else:
            self.tn = int(np.floor(crop/100*dwi.shape[-1]))

        if kernel is None:
            nvols = dwi.shape[-1]
            p = np.arange(3, nvols, 2)
            pf = np.where(p**3 >= nvols)[0][0]
            defaultKernel = p[pf]
            kernel = np.array([defaultKernel,defaultKernel,defaultKernel])
        else:
            kernel = np.array(kernel)

        kernel = kernel+np.mod(kernel,2)-1
        self.kernel = kernel.astype(int)
        k = self.kernel // 2
        
        if patchtype is None:
            self.patchtype = 'box'
        else:
            self.patchtype = patchtype

        if patchsize is None:
            self.patchsize = np.prod(kernel)
        else:
            self.patchsize = patchsize

        if self.patchtype == 'box':
            self.pos_distances = None
        elif self.patchtype == 'nonlocal':
            if self.patchsize >= np.prod(kernel):
                print('Warning: selecting sane default nonlocal patch size')
                self.patchsize = np.floor(np.prod(kernel) - 0.2*np.prod(kernel))
                if self.patchsize <= dwi.shape[-1]:
                    self.patchsize = dwi.shape[-1] + 1
                self.patchsize = int(self.patchsize)
            pi, pj, pk = np.where(np.ones((self.kernel)))
            patchcoords = np.vstack((pi,pj,pk)).T
            self.pos_distances = np.sum((patchcoords - k)**2, axis=1)
                    
        if self.coil:
            pwidth = (k[0],k[0]), (k[1],k[1]), (k[2],k[2]), (0,0), (0,0)
        else:
            pwidth = (k[0],k[0]), (k[1],k[1]), (k[2],k[2]), (0,0)
        self.pwidth = pwidth
        self.dwi = np.pad(dwi, pad_width=pwidth, mode='wrap')

        if pmap is None:
            self.pmap = None
        else:
            self.pmap = np.pad(pmap, pad_width=pwidth[:][:-2],mode='wrap')
            
        print('Denoising data with parameters:')
        print('kernel     = ' + str(self.kernel))
        print('patch type = ' + str(self.patchtype))
        print('patch size = ' + str(self.patchsize))
        print('shrinkage  = ' + str(self.shrink))
        print('algorithm  = ' + str(self.algo))
        print('coil level = ' + str(self.coil))
        print('rotation   = ' + str(self.rotation))
        print('cropping   = ' + str(self.tn) + ' (' + str(crop) + '%)')
            
    def box_patch(self, dwi, coords):
        # extracts a patch of size kx x ky x kz from the padded input at specified coords
        k = self.kernel // 2
        X = dwi[coords[0]-k[0]:coords[0]+k[0]+1, coords[1]-k[1]:coords[1]+k[1]+1, coords[2]-k[2]:coords[2]+k[2]+1, ...]
        return X

    def normalize(self, im):
        im = ((im - np.min(im)) * (1/(np.max(im) - np.min(im)) * 1.0))
        return im

    def refine(self, X_tmp, coords):
        refval = X_tmp[np.prod(self.kernel)//2,...]
        if self.coil:
            int_distances = 1/(X_tmp.shape[-1]*X_tmp.shape[-2]) * np.sum((X_tmp - refval)**2, axis=(1,2))
        else:
            int_distances = 1/(X_tmp.shape[-1]) * np.sum((X_tmp - refval)**2, axis=1)

        wdists = (self.pos_distances * int_distances)
        iidx = np.argpartition(wdists, self.patchsize)[:self.patchsize]
        minind = np.argmin(self.pos_distances[iidx])
        return iidx, minind

    def unpad(self, x, pad_width):
        slices = []
        for c in pad_width:
            e = None if c[1] == 0 else -c[1]
            slices.append(slice(c[0], e))
        return x[tuple(slices)]

    def padded_sampling(self, mask):
        # outputs a grid x, y, z of which coordinates to loop over when processing
        sx, sy, sz = mask.shape
        k = self.kernel // 2
        mask[:k[0],:,:] = 0
        mask[sx-k[0]:,:,:] = 0
        mask[:,:k[1],:] = 0
        mask[:,sy-k[1]:,:] = 0
        mask[:,:,:k[2]] = 0
        mask[:,:,sz-k[2]:] = 0
        self.mask = mask
        x, y, z = np.where(mask==1)
        return x.astype(int), y.astype(int), z.astype(int)

    def sampling(self, mask):
        # outputs a grid x, y, z of which coordinates to loop over when processing
        x, y, z = np.where(mask==1)
        return x.astype(int), y.astype(int), z.astype(int)

    def eig_shrink(self, vals, gamma):
        t = 1 + np.sqrt(gamma)
        s = np.zeros((vals.shape))
        x = vals[vals > t]
        s[vals > t] = np.sqrt((x**2 - gamma - 1)**2 - 4*gamma) / x
        return np.diag(s)

    def jsvd(self, X):
        u,vals,v = jsp.linalg.svd(X, full_matrices=False, compute_uv=True, check_finite=False, overwrite_a=True)
        return u,vals,v
    
    def denoise(self, coords):
        Xbox = self.box_patch(self.dwi, coords)

        if self.coil:
            X = np.reshape(Xbox,(np.prod(self.kernel), Xbox.shape[-2], Xbox.shape[-1]), order='F')
        else:
            X = np.reshape(Xbox,(np.prod(self.kernel), Xbox.shape[-1]), order='F')

        if self.patchtype == 'nonlocal':
            Xn = self.normalize(X)
            nonlocalinds, minind = self.refine(Xn, coords)
            X = X[nonlocalinds,...]
            if self.pmap is not None:
                P = self.pmap[coords[0],coords[1],coords[2]].astype(int)
                X = X[:,np.arange(P).astype(int),...]
        else:
            minind = np.prod(self.kernel)//2

        if self.coil:
            X = X.reshape(self.patchsize*X.shape[-2], X.shape[-1], order='F')
        
        M = X.shape[0]
        N = X.shape[1]
        Mp = np.min((M,N))
        Np = np.max((M,N))
        
        if M < N:  
            X = np.conj(X).T
        
        #try:
        #u,vals,v = jit(self.jsvd)(X)
        u,vals,v = scipy.linalg.svd(X, full_matrices=False, overwrite_a=True, check_finite=False)
        irot = None
        
        vals = (vals**2).astype('float32')
        order = np.argsort(vals)[::-1]
        u = u[:,order]
        v = v[:,order]
        vals = vals[order]

        ptn = np.arange(0,Mp-self.tn).T
        p = np.arange(0,Mp).T

        csum = np.cumsum(vals[::-1])[::-1]
        if self.algo == 'veraart':
            sigmasq_1 = csum / ((Mp-p)*Np)
            rangeMP = 4*np.sqrt((Mp-ptn)*(Np-self.tn))
        elif self.algo == 'cordero-grande':
            sigmasq_1 = csum / ((Mp-p)*(Np-p))
            rangeMP = 4*np.sqrt((Mp-ptn)*(Np-ptn))
        elif self.algo == 'jespersen':
            sigmasq_1 = csum / ((Mp-p)*(Np-p))
            rangeMP = 4*np.sqrt((Np-self.tn)*(Mp))

        rangeData = vals[:Mp-self.tn] - vals[Mp-self.tn-1]
        sigmasq_2 = rangeData / rangeMP

        t = np.where(sigmasq_2 < sigmasq_1[:Mp-self.tn])
        if t[0].any():
            t = t[0][0]
        else:
            t = Mp-1
        
        sigma = np.sqrt(sigmasq_1[t])
        npars = t
        if self.shrink == 'threshold':
            vals_t = vals.copy()
            vals_t[t:] = 0
            s = np.array(u) @ np.diag(np.sqrt(vals_t)) @ np.array(v)
        else:
            vals_norm = np.sqrt(vals)/(np.sqrt(Mp)*sigma)
            vals_frobnorm = (np.sqrt(Mp)*sigma) * self.eig_shrink(vals_norm, Np/Mp)
            s = np.array(u) @ vals_frobnorm @ np.array(v) 
        
        if M<N:
            s = np.conj(s).T

        if self.coil:
            signal = np.squeeze(s[minind::self.patchsize,:])
            if self.pmap is not None:
                signal = np.vstack((signal, np.zeros((self.dwi.shape[-2]-P,self.dwi.shape[-1]))))
        else:    
            signal = np.squeeze(s[minind, :]) 

        return signal, sigma, npars, irot, v

    def denoise_g(self, coords):
        Xbox = self.box_patch(self.dwi, coords)
        X = Xbox.reshape(np.prod(self.kernel),Xbox.shape[-2],Xbox.shape[-1],order='F').transpose(0,2,1).reshape(np.prod(self.kernel)*Xbox.shape[-1],Xbox.shape[-2],order='F')  

        M = X.shape[0]
        N = X.shape[1]
        Mp = np.min((M,N))
        Np = np.max((M,N))
        
        if M > N:  
            X = np.conj(X).T

        u,vals,v = scipy.linalg.svd(X, full_matrices=False, overwrite_a=True, check_finite=False)
        vs = np.diag(vals) @ np.array(v)
        vals = (vals**2)
        
        # order = np.argsort(vals)[::-1]
        # u = u[:,order]
        # v = v[:,order]
        # vals = vals[order]

        tn=5   
        ptn = np.arange(0,Mp-tn).T
        p = np.arange(0,Mp).T

        csum = np.cumsum(vals[::-1])[::-1]
        if self.algo == 'veraart':
            sigmasq_1 = csum / ((Mp-p)*Np)
            rangeMP = 4*np.sqrt((Mp-ptn)*(Np-tn))
        elif self.algo == 'cordero-grande':
            sigmasq_1 = csum / ((Mp-p)*(Np-p))
            rangeMP = 4*np.sqrt((Mp-ptn)*(Np-ptn))
        elif self.algo == 'jespersen':
            sigmasq_1 = csum / ((Mp-p)*(Np-p))
            rangeMP = 4*np.sqrt((Np-tn)*(Mp))

        rangeData = vals[:Mp-tn] - vals[Mp-tn-1]
        sigmasq_2 = rangeData / rangeMP

        t = np.where(sigmasq_2 < sigmasq_1[:Mp-tn])
        if t[0].any():
            t = t[0][0]
            pass

        sigma = np.sqrt(sigmasq_1[t])
        npars = t
        
        if self.shrink == 'threshold':
            vals_t = vals.copy()
            vals_t[t:] = 0
            s = np.array(v).T @ np.diag(np.sqrt(vals_t))
        else:
            vals_norm = np.sqrt(vals)/(np.sqrt(Mp)*sigma)
            vals_frobnorm = (np.sqrt(Mp)*sigma) * self.eig_shrink(vals_norm, Np/Mp)
            s = np.array(v).T @ vals_frobnorm
        
        if M>N:
            s = np.conj(s).T

        imd = vs.reshape(Xbox.shape[-2],self.kernel[0],self.kernel[1],self.kernel[2], Xbox.shape[-1], order='F').transpose(1,2,3,0,4)
        irot = imd[self.kernel[0]//2,self.kernel[1]//2,self.kernel[2]//2,:,:]

        s = s.reshape(Xbox.shape[-2],self.kernel[0],self.kernel[1],self.kernel[2], Xbox.shape[-1], order='F').transpose(1,2,3,0,4)
        signal = s[self.kernel[0]//2,self.kernel[1]//2,self.kernel[2]//2,:,:]
        
        return signal, sigma, npars, irot, np.array(u)

    def unrotate(self,coords, VV):
        Xbox = self.box_patch(self.dwi, coords)
        X = Xbox.reshape(np.prod(self.kernel),Xbox.shape[-2],Xbox.shape[-1],order='F').transpose(0,2,1).reshape(np.prod(self.kernel)*Xbox.shape[-1],Xbox.shape[-2],order='F')  
        if X.shape[1] == self.dwi.shape[-2]:
            X = np.conj(X).T
        datap = VV @ X
        datap = datap.reshape(self.dwi.shape[-2], self.kernel[0], self.kernel[1], self.kernel[2], self.dwi.shape[-1],order='F').transpose(1,2,3,0,4)
        datap = datap[self.kernel[0]//2,self.kernel[1]//2,self.kernel[2]//2,:,:]
        return datap

    def undo_rotation(self, VV):
        if self.coil:
            sx, sy, sz, sc, N = self.dwi.shape
        else:
            sx, sy, sz, N = self.dwi.shape
        mask = np.ones((sx,sy,sz))

        x, y, z = self.padded_sampling(mask)
        xsize = int(x.size)
        coords = np.vstack((x,y,z))
        
        inputs = tqdm(range(0, xsize))
        
        VV_pad = np.pad(VV,pad_width=self.pwidth,mode='reflect')

        signal = Parallel(n_jobs=-1, prefer='processes')\
            (delayed(self.unrotate)(coords[:,i],VV_pad[coords[0,i],coords[1,i],coords[2,i],:,:]) for i in inputs)

        if self.coil:
            Signal = np.zeros((sx, sy, sz, sc, N)).astype(complex)
        else:
            Signal = np.zeros((sx, sy, sz, N)).astype(self.dwi.dtype)

        for nn in inputs:
            Signal[x[nn], y[nn], z[nn], ...] = signal[nn]
        Signal = self.unpad(Signal, self.pwidth)
        return Signal


    def process(self):
        if self.coil:
            sx, sy, sz, sc, N = self.dwi.shape
        else:
            sx, sy, sz, N = self.dwi.shape
        mask = np.ones((sx,sy,sz))

        x, y, z = self.padded_sampling(mask)
        xsize = int(x.size)
        coords = np.vstack((x,y,z))
        
        inputs = tqdm(range(0, xsize))
        num_cores = multiprocessing.cpu_count()
        
        # # # parallel
        if self.rotation:
            signal, sigma, npars, irot, v = zip(*Parallel(n_jobs=-1, prefer='processes')\
                (delayed(self.denoise_g)(coords[:,i]) for i in inputs))
        else:
            signal, sigma, npars, irot, v = zip(*Parallel(n_jobs=-1, prefer='processes')\
                (delayed(self.denoise)(coords[:,i]) for i in inputs))

      #  serial
        # k = self.kernel // 2
        # # for t in inputs:
        # crds = np.array([70+k[0], 35+k[1], 0+k[2]])
        # a, b, c, d, e = self.denoise_g(crds)
        # print(c)
        # # print(np.shape(e))
        # # import pdb
        # # pdb.set_trace()
           
        # reconstruct original data matrix
        Sigma = np.zeros((sx, sy, sz))
        Npars = np.zeros((sx, sy, sz))
        if self.coil:
            Signal = np.zeros((sx, sy, sz, sc, N)).astype(complex)
        else:
            Signal = np.zeros((sx, sy, sz, N)).astype(self.dwi.dtype)
        if self.rotation:
            Rotv = np.zeros((sx, sy, sz, sc, sc)).astype(complex)
            Irot = np.zeros((sx, sy, sz, sc, N)).astype(complex)

        for nn in inputs:
            Sigma[x[nn], y[nn], z[nn]] = sigma[nn]
            Npars[x[nn], y[nn], z[nn]] = npars[nn]
            Signal[x[nn], y[nn], z[nn], ...] = signal[nn]
            if self.rotation:
                Rotv[x[nn], y[nn], z[nn], ...] = v[nn]
                Irot[x[nn], y[nn], z[nn], ...] = irot[nn]
                
        Signal = self.unpad(Signal, self.pwidth)
        if self.coil:
            Npars = self.unpad(Npars, self.pwidth[:][:-2])
            Sigma = self.unpad(Sigma, self.pwidth[:][:-2])
        else:
            Npars = self.unpad(Npars, self.pwidth[:][:-1])
            Sigma = self.unpad(Sigma, self.pwidth[:][:-1])
        if self.rotation:
            Rotv = self.unpad(Rotv, self.pwidth)
            Irot = self.unpad(Irot, self.pwidth)
        else:
            Rotv = None
            Irot = None

        return Signal, Sigma, Npars, Irot, Rotv

def denoise(img, kernel=None, patchtype=None, patchsize=None, shrinkage=None, algorithm=None, rotation=None, pmap=None, crop=None):
    mp = MP(img, kernel, patchtype, patchsize, shrinkage, algorithm, rotation, pmap, crop)
    Signal, Sigma, Npars, Irot, Rotu = mp.process()
    return Signal, Sigma, Npars, Irot, Rotu

def unrotate(img, VV, kernel=None,patchtype=None):
    mp = MP(img, kernel, patchtype, patchsize=None, shrinkage=None, algorithm=None, rotation=None, pmap=None, crop=None)
    Signal = mp.undo_rotation(VV)
    return Signal


# if __name__ == "__main__":
#     denoise()
    
