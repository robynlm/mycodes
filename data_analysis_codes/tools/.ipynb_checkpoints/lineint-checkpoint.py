"""This module provides the LineIntegrate class.

Copyright (C) 2024  Robyn L. Munoz

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

You may contact the author at : robyn.munoz@yahoo.fr
"""

import numpy as np

class LineIntegrate:
    """Class to compute the numerical distance of a line in a grid.
    
    This class provides for each grid position traversed by the line 
    the size of the segment contained within.
    
    """
    def __init__(self, xdomain, ydomain, zdomain, idx_start):
        """Define grid and sphere center
        
        Parameters :
            xdomain : (3) list of float, [xmin, xmax, dx]
            ydomain : (3) list of float, [ymin, ymax, dy]
            zdomain : (3) list of float, [zmin, zmax, dz]
            idx_start : (3) list of int, indices of strating location
                        (from where the lines start)
        """
        
        # grid
        self.xmin, self.xmax, self.dx = xdomain
        self.ymin, self.ymax, self.dy = ydomain
        self.zmin, self.zmax, self.dz = zdomain
        self.Lx = self.xmax - self.xmin
        self.Ly = self.ymax - self.ymin
        self.Lz = self.zmax - self.zmin
        self.Nx = self.Lx/self.dx
        self.Ny = self.Ly/self.dy
        self.Nz = self.Lz/self.dz
        self.d1x = np.arange(self.xmin, self.xmax, self.dx)
        self.d1y = np.arange(self.ymin, self.ymax, self.dy)
        self.d1z = np.arange(self.zmin, self.zmax, self.dz)
        self.cell_diag = np.sqrt(self.dx**2 + self.dy**2 + self.dz**2)
        
        # starting coordinates and indices
        self.ixC, self.iyC, self.izC = idx_start
        self.xC = self.d1x[self.ixC]
        self.yC = self.d1y[self.iyC]
        self.zC = self.d1z[self.izC]        
        
    def safe_division(self, numerator, denominator):
        """
        """
        num0 = numerator == 0
        den0 = denominator == 0
        if isinstance(numerator, float):
            ratio = 0.0 if den0 else numerator/denominator
        else: # it's an array
            ratio = np.divide(numerator, denominator, 
                              out = np.zeros_like(numerator), 
                              where = ~den0)
            ratio[np.where(np.logical_and(num0, den0))] = 1.0
        return ratio 
        
    def radius3D(self, x, y, z):
        """Compute 3D radius.
        
        Parameters :
            x : float, location along x
            y : float, location along y
            z : float, location along z
            
        Returns :
            radius : float
        """
        return np.sqrt((x)**2 + (y)**2 + (z)**2)
        
    def radius2D(self, x, y):
        """Compute 2D radius.
        
        Parameters :
            x : float, location along x
            y : float, location along y
            
        Returns :
            radius : float
        """
        return np.sqrt((x)**2 + (y)**2)
    
    def check(self, i, theta, phi, A, B):
        theta_deg = theta*180/np.pi
        phi_deg = phi*180/np.pi
        if abs(A/B - 1)>1e-8:
            print('WARNING:' + i 
                  + ' theta={:.2f}, phi={:.2f},'.format(theta_deg, phi_deg)
                  + ' {:.5f}!={:.5f}'.format(A, B))
        
    def segments_idx_end(self, idx_end):
        """
        
        Parameters :
            idx_end : (3) list of int, indices of 
                      end point of the line
                      
        Returns :
            indices : (?, 3) list of tuple of int
                      ? is the number of grid points crossed by the line
                      the tuple of 3 ints is the grid points' indices
            segments : (?) np.array of float
                       ? is the number of grid points crossed by the line
                       the float is the size of corresponding segment
            theta : float, inclination angle of line
            phi : float, azimuth angle of line
        """
        if idx_end == (self.ixC, self.iyC, self.izC):
            # if starting point = end point
            return [idx_end], [0.0], 0.0, 0.0
        else:
            # === Cartesian coordinates
            # Shift so that starting location is in the center
            x = -self.xC + self.d1x[idx_end[0]]
            y = -self.yC + self.d1y[idx_end[1]]
            z = -self.zC + self.d1z[idx_end[2]]
            if x >= self.xmax:
                x = x - self.Lx
            if y >= self.ymax:
                y = y - self.Ly
            if z >= self.zmax:
                z = z - self.Lz

            # === Spherical coordinates
            # radius
            r = self.radius3D(x, y, z)
            # inclination
            theta = np.arccos(self.safe_division(z, r))
            if r==0:
                theta = 0
            # azimuth
            phi = np.arccos(self.safe_division(x, self.radius2D(x, y)))
            if y<0:
                phi = 2*np.pi - phi
            if x==0 and y==0:
                phi = 0

            # === Check coord location
            xc = r * np.cos(phi) * np.sin(theta)
            yc = r * np.sin(phi) * np.sin(theta)
            zc = r * np.cos(theta)
            self.check('x', theta, phi, x, xc)
            self.check('y', theta, phi, y, yc)
            self.check('z', theta, phi, z, zc)

            # === Max indices
            buffer = 4
            rmax = r + self.cell_diag
            nr = r + buffer * self.cell_diag
            nnr = nr + buffer * self.cell_diag
            rx, ry, rz = [i for _,i in sorted(zip(np.array([x,y,z]), 
                                                  np.array([r,nr,nnr])))]
            xmax = rx * np.cos(phi) * np.sin(theta)
            ymax = ry * np.sin(phi) * np.sin(theta)
            zmax = ry * np.cos(theta)
            imax = [int(abs(i/self.dx))+buffer,
                    int(abs(i/self.dy))+buffer,
                    int(abs(i/self.dz))+buffer]

            # === Main calculation
            indices, segments = self.idx_along_r(theta, phi, imax, rmax, 
                                                 idx_end, x, y, z)
            # shift around back to starting location
            findices = []
            for idx in indices:
                ix, iy, iz = idx
                ix = ix - self.ixC
                iy = iy - self.iyC
                iz = iz - self.izC
                if ix < 0:
                    ix = self.Nx + ix
                if iy < 0:
                    iy = self.Ny + iy
                if iz < 0:
                    iz = self.Nz + iz
                if ix > self.Nx - 1:
                    ix = ix - self.Nx
                if iy > self.Ny - 1:
                    iy = iy - self.Ny
                if iz > self.Nz - 1:
                    iz = iz - self.Nz
                findices += [(ix, iy, iz)]

            # cut data to end data point
            try:
                i = findices.index(idx_end) + 1
            except:
                print(x, y, z)
                print(indices)
                print(findices)
                i = findices.index(idx_end) + 1
            final_indices = findices[:i]
            final_segments = list(segments[:i])
            
            # fix last weight
            numerical_radius = np.sum(final_segments)
            if r < numerical_radius:
                print('need to fix last weight')
                final_segments += [r - radius[i]]
                final_indices += [findices[i - 1]]
            elif r > numerical_radius:
                print('need to fix last weight')
                final_segments += [r - radius[i]]
                final_indices += [findices[i]]

            return final_indices, final_segments, theta, phi
        
    def idx_along_r(self, theta, phi, imax, rmax,
                    idx_end, xend, yend, zend):
        """
        * : grid indice
        -|: grid lines
        
         * | * | *
        ---|---|---
         * | * | *
        ---|---|---
         * | * | *
        
        """
        
        #================================================================
        # Get the coordinates where the line intersects with the grid
        #================================================================
        
        # angles
        signx = -1 if phi > np.pi/2 and phi<3*np.pi/2 else 1
        signy =  1 if phi <= np.pi else -1
        signz =  1 if theta <= np.pi /2 else -1
        good_xangle = theta!=0 and theta!=np.pi and phi!=np.pi/2 and phi!=3*np.pi/2
        good_yangle = theta!=0 and theta!=np.pi and phi!=0 and phi!=np.pi
        good_zangle = theta!=np.pi/2

        # intersect with x lines
        dirx_coordx = np.array([signx * self.dx * (0.5 + i)
                                     for i in range(imax[0]) if good_xangle])
        dirx_coordy = dirx_coordx * np.tan(phi)
        dirx_coordz = self.safe_division(self.radius2D(dirx_coordx,
                                                       dirx_coordy),
                                         np.tan(theta))

        # intersect with y lines
        diry_coordy = np.array([signy * self.dy * (0.5 + i)
                                     for i in range(imax[1]) if good_yangle])
        diry_coordx = self.safe_division(diry_coordy, np.tan(phi))
        diry_coordz = self.safe_division(self.radius2D(diry_coordx,
                                                       diry_coordy),
                                         np.tan(theta))

        # intersect with z lines
        dirz_coordz = np.array([signz * self.dz * (0.5 + i)
                                     for i in range(imax[2]) if good_zangle])
        dirz_coordx = dirz_coordz * np.cos(phi) * np.tan(theta)
        dirz_coordy = dirz_coordz * np.sin(phi) * np.tan(theta)

        #================================================================
        # Calculate size of segments between intersections
        #================================================================
    
        # put all the coordinates together
        coordx = [0.0] + list(dirx_coordx) + list(diry_coordx) + list(dirz_coordx)
        coordy = [0.0] + list(dirx_coordy) + list(diry_coordy) + list(dirz_coordy)
        coordz = [0.0] + list(dirx_coordz) + list(diry_coordz) + list(dirz_coordz)
        
        # calc radius of each coordinate
        radius_dirx = [self.radius3D(xi, yi, zi) 
                       for xi, yi, zi in zip(dirx_coordx, dirx_coordy, dirx_coordz)]
        radius_diry = [self.radius3D(xi, yi, zi) 
                       for xi, yi, zi in zip(diry_coordx, diry_coordy, diry_coordz)]
        radius_dirz = [self.radius3D(xi, yi, zi) 
                       for xi, yi, zi in zip(dirz_coordx, dirz_coordy, dirz_coordz)]
        radius = [0.0] + radius_dirx + radius_diry + radius_dirz
        
        # sort the points according to radius size
        coordx = [i for _,i in sorted(zip(radius, coordx))]
        coordy = [i for _,i in sorted(zip(radius, coordy))]
        coordz = [i for _,i in sorted(zip(radius, coordz))]
        radius.sort()
        
        # set radius cutoff as min(max(radius in different directions))
        radius_cutoff = []
        for radius_dir in [radius_dirx, radius_diry, radius_dirz]:
            if radius_dir != []:
                radius_cutoff += [np.max(radius_dir)]
        radius_cutoff = np.min(radius_cutoff)
        if radius_cutoff < rmax:
            radius_cutoff = rmax
        
        # remove repeated radius values and exclude points > cutoff
        coordx_retained, coordy_retained, coordz_retained = [], [], []
        radius_retained = []
        for i, r in enumerate(radius):
            # TO DO: I don't like using round, change that
            if np.round(r, 3) not in radius_retained and np.round(r, 3) < radius_cutoff:
                radius_retained += [np.round(r, 3)]
                coordx_retained += [coordx[i]]
                coordy_retained += [coordy[i]]
                coordz_retained += [coordz[i]]
        coordx = coordx_retained
        coordy = coordy_retained
        coordz = coordz_retained
        coord  = [coordx, coordy, coordz]
        radius = radius_retained

        # calculate segments
        segments = np.array([np.sqrt((coordx[i+1] - coordx[i])**2
                                     + (coordy[i+1] - coordy[i])**2
                                     + (coordz[i+1] - coordz[i])**2) 
                             for i in range(len(coordx)-1)])
        
        # check
        numerical_radius = np.sum(segments)
        self.check('1', theta, phi, numerical_radius, radius[-1])

        #================================================================
        # Get associated index in the grid
        #================================================================

        # the average between the intersections 
        # will be closer to the position of the data points
        # x indice
        av_coordx = [np.average([coordx[i], coordx[i+1]]) 
                     for i in range(len(coordx) - 1)]
        xidx = [np.argmin(abs(self.d1x - i)) for i in av_coordx] 
        
        # y indice
        av_coordy = [np.average([coordy[i], coordy[i+1]]) 
                     for i in range(len(coordy) - 1)]
        yidx = [np.argmin(abs(self.d1y - i)) for i in av_coordy]
        
        # z indice
        av_coordz = [np.average([coordz[i], coordz[i+1]]) 
                     for i in range(len(coordz) - 1)]
        zidx = [np.argmin(abs(self.d1z - i)) for i in av_coordz]
        
        # put the indices together
        indices = [(ix, iy, iz) for ix, iy, iz in zip(xidx, yidx, zidx)]
        
        return indices, segments