"""
Author: Robyn Munoz
Date: 4th Jan 2021

This class applies a 4th order finite difference approximation of a spatial derivative.
The formulas used here are from A.Yew (2011), and assume the boundaries to be periodic.
"""

import numpy as np

class FD4():
    def backward(self, f, i, Delta):
        return ((25/12)*f[i] 
                + (-4)*f[i-1] 
                + (3)*f[i-2] 
                + (-4/3)*f[i-3] 
                + (1/4)*f[i-4])/Delta
    def centered(self, f, i, Delta):
        return ((1/12)*f[i-2] 
                + (-2/3)*f[i-1] 
                + (2/3)*f[i+1] 
                + (-1/12)*f[i+2])/Delta
    def forward(self, f, i, Delta):
        return ((-25/12)*f[i] 
                + (4)*f[i+1] 
                + (-3)*f[i+2] 
                + (4/3)*f[i+3] 
                + (-1/4)*f[i+4])/Delta
    
class FD6():
    def backward(self, f, i, Delta):
        return ((49/20)*f[i] 
                + (-6)*f[i-1] 
                + (15/2)*f[i-2] 
                + (-20/3)*f[i-3] 
                + (15/4)*f[i-4] 
                + (-6/5)*f[i-5] 
                + (1/6)*f[i-6])/Delta
    def centered(self, f, i, Delta):
        return ((-1/60)*f[i-3] 
                + (3/20)*f[i-2] 
                + (-3/4)*f[i-1] 
                + (3/4)*f[i+1] 
                + (-3/20)*f[i+2] 
                + (1/60)*f[i+3])/Delta
    def forward(self, f, i, Delta):
        return ((-49/20)*f[i] 
                + (6)*f[i+1] 
                + (-15/2)*f[i+2] 
                + (20/3)*f[i+3] 
                + (-15/4)*f[i+4] 
                + (6/5)*f[i+5] 
                + (-1/6)*f[i+6])/Delta   

class FD_Class():
    def __init__(self, 
                 dx, dy, dz, 
                 periodic_boundary=True, 
                 order6=False):
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.pBound = periodic_boundary
        
        if order6:
            self.FD = FD6()
            self.mask = 3
        else:
            self.FD = FD4()
            self.mask = 2
            
    def Dt(self, f, dt):
        return np.array(
            [self.FD.backward(f, it, dt) 
             for it in range(self.mask*2,len(f))])
    
    def D3x(self, f): # 3D Derivative along x
        if self.pBound:
            flong = np.concatenate(
                (f[-self.mask:, :, :], f, 
                 f[:self.mask, :, :]), axis=0)
            return np.array(
                [self.FD.centered(flong, ix, self.dx) 
                 for ix in range(
                     self.mask, len(f)+self.mask)])
        else:
            N = len(f)
            LHS = np.array(
                [self.FD.forward(f, ix, self.dx) 
                 for ix in range(0, self.mask)])
            central_part = np.array(
                [self.FD.centered(f, ix, self.dx) 
                 for ix in range(self.mask, N-self.mask)])
            RHS = np.array(
                [self.FD.backward(f, ix, self.dx) 
                 for ix in range(N-self.mask, N)])
            return np.concatenate(
                (LHS, central_part, RHS), axis=0)
    
    def D3y(self, f): # 3D Derivative along y
        N = len(f)
        if self.pBound:
            flong = np.concatenate(
                (f[:, -self.mask:, :], f,
                 f[:, :self.mask, :]), axis=1)
            return np.array(
                [[self.FD.centered(
                    flong[ix,:,:], iy, self.dy) 
                  for iy in range(self.mask, N+self.mask)] 
                 for ix in range(N)])
        else:
            LHS = np.array(
                [[self.FD.forward(f[ix,:,:], iy, self.dy) 
                  for iy in range(0, self.mask)] 
                 for ix in range(N)])
            central_part = np.array(
                [[self.FD.centered(f[ix,:,:], iy, self.dy) 
                  for iy in range(self.mask, N-self.mask)] 
                 for ix in range(N)])
            RHS = np.array(
                [[self.FD.backward(f[ix,:,:], iy, self.dy) 
                  for iy in range(N-self.mask, N)] 
                 for ix in range(N)])
            return np.concatenate(
                (LHS, central_part, RHS), axis=1)
    
    def D3z(self, f): # 3D Derivative along z
        N = len(f)
        if self.pBound:
            flong = np.concatenate(
                (f[:, :, -self.mask:], f, 
                 f[:, :, :self.mask]), axis=2)
            return np.array(
                [[[self.FD.centered(
                    flong[ix,iy,:], iz, self.dz) 
                   for iz in range(self.mask, N+self.mask)] 
                  for iy in range(N)] 
                 for ix in range(N)])
        else:
            LHS = np.array(
                [[[self.FD.forward(f[ix,iy,:], iz, self.dz) 
                   for iz in range(0, self.mask)] 
                  for iy in range(N)] 
                 for ix in range(N)])
            central_part = np.array(
                [[[self.FD.centered(f[ix,iy,:], iz, self.dz) 
                   for iz in range(self.mask, N-self.mask)] 
                  for iy in range(N)] 
                 for ix in range(N)])
            RHS = np.array(
                [[[self.FD.backward(f[ix,iy,:], iz, self.dz) 
                   for iz in range(N-self.mask, N)] 
                  for iy in range(N)] 
                 for ix in range(N)])
            return np.concatenate(
                (LHS, central_part, RHS), axis=2)
    
    
    # The following functions apply the 3 derivatives to a scalar or a tensor
    def D3_scalar(self, f): # \partial_i (f)
        return np.array(
            [self.D3x(f), self.D3y(f), self.D3z(f)])
    def D3_tensor1(self, f): # \partial_i (f_{j}) or \partial_i (f^{j})
        dim = len(f)
        return np.array(
            [[self.D3x(f[j]) for j in range(dim)],
             [self.D3y(f[j]) for j in range(dim)],
             [self.D3z(f[j]) for j in range(dim)]])
    def D3_tensor2(self, f): # \partial_i (f_{kj}) or \partial_i (f^{kj})
        dim = len(f)
        return np.array(
            [[[self.D3x(f[k, j]) 
               for j in range(dim)] 
              for k in range(dim)],
             [[self.D3y(f[k, j]) 
               for j in range(dim)] 
              for k in range(dim)],
             [[self.D3z(f[k, j]) 
               for j in range(dim)] 
              for k in range(dim)]])
