import numpy as np
from . import ReadingTools as RRead

class TA_Class:
    def __init__(self, param, Lin, OD_Loc):
    
        self.param = param
        self.dx = param['dx']
        self.L = param['Lx']
        self.N = param['Nx']
        self.Lin = Lin
        
        self.P4 = np.pi/4
        self.P2 = np.pi/2
        
        self.OD_loc = OD_Loc
        self.ixOD = OD_Loc[0][0]
        self.iyOD = OD_Loc[1][0]
        self.izOD = OD_Loc[2][0]
        self.xOD = self.Lin.d1x[self.ixOD]
        self.yOD = self.Lin.d1y[self.iyOD]
        self.zOD = self.Lin.d1z[self.izOD]
        self.xyzOD = np.array([self.xOD, self.yOD, self.zOD])
        self.radius_in_grid = np.sqrt((self.Lin.d3x-self.xOD)**2
                                      + (self.Lin.d3y-self.yOD)**2
                                      + (self.Lin.d3z-self.zOD)**2)
        
        r_possible = []
        for ix in self.Lin.d1x:
            for iy in self.Lin.d1y:
                for iz in self.Lin.d1z:
                    rnew = np.sqrt(ix**2+iy**2+iz**2)
                    if rnew < np.sqrt(3)*self.L:
                        r_possible += [rnew]
        r_possible = np.array(list(dict.fromkeys(r_possible)))
        r_possible = np.append(-r_possible, r_possible)
        r_possible.sort()
        self.r_possible = r_possible
        
    def xyz_periodic(self, i_all, cut=True, boundary_crossed=False):
        if isinstance(i_all, (int, float, np.generic)):
            i_all = [i_all]
        xyz_return = []
        for i in i_all:
            if i<self.N:
                if cut and self.Lin.d1x[i]>self.L/4:
                    xyz_return += [self.Lin.d1x[0] 
                                   - (self.L/2-self.Lin.d1x[i])]
                    boundary_crossed=True  
                    #I assume that the idx are sorted according to radius
                elif boundary_crossed:
                    xyz_return += [self.Lin.d1x[0] 
                                   - (self.L/2-self.Lin.d1x[i])]
                else:
                    xyz_return += [self.Lin.d1x[i]]
            else:
                xyz_return += [self.Lin.d1x[0] - ((self.N-1)-i)*self.dx]
                
        if len(xyz_return)==1:
            return xyz_return[0], boundary_crossed
        else:
            return np.array(xyz_return), boundary_crossed
        
    def radius_idx(self, idx, cut=True):
        x = np.array([self.Lin.d1x[i[0]] for i in idx])
        y = np.array([self.Lin.d1x[i[1]] for i in idx])
        z = np.array([self.Lin.d1x[i[2]] for i in idx])
        return self.radius(x, y, z)
    
    def radiusOD(self, x, y, z):
        r = np.sqrt((x-self.xOD)**2 + (y-self.yOD)**2 + (z-self.zOD)**2)
        return r
    
    def radius(self, x, y, z):
        r = np.sqrt((x)**2 + (y)**2 + (z)**2)
        return r
    
    def radius2DOD(self, x, y):
        r = np.sqrt((x-self.xOD)**2 + (y-self.yOD)**2)
        return r
    
    def radius2D(self, x, y):
        r = np.sqrt((x)**2 + (y)**2)
        return r
        
    def find2min(self, V):
        iv1 = np.argmin(V)
        V2 = np.append(V[:iv1],V[iv1+1:])
        iv2 = np.argmin(V2)
        if iv2>=iv1:
            iv2 += 1
        if iv2==iv1-1:
            iv2 = iv1
        return np.array([iv1, iv2])
    
    def find4min(self, V):
        idx_min = []
        for i in range(2):
            iv = np.argmin(V)
            V[iv] += 1e5
            if iv<len(V)-1:
                ivp = iv-1+2*np.argmin([V[iv-1], V[iv+1]])
            else:
                ivp = iv-1 if V[iv-1]<V[len(V)-iv] else len(V)-iv
            V[ivp] += 1e5
            idx_min += list(np.sort([iv,ivp]))
        return idx_min
    
    
    def closest_idx(self, coord):
        return np.argmin(abs(self.Lin.d1x - coord))

    def idx_along_r(self, theta, phi, imax, rmax,
                    idx_grid_point_considered, xwanted, ywanted, zwanted):
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
        # Get the coordinates where the  line intersects with the grid
        #================================================================
        
        # angles
        signx = -1 if phi > np.pi/2 and phi<3*np.pi/2 else 1
        signy =  1 if phi <= np.pi else -1
        signz =  1 if theta <= np.pi /2 else -1
        good_xangle = theta!=0 and theta!=np.pi and phi!=np.pi/2 and phi!=3*np.pi/2
        good_yangle = theta!=0 and theta!=np.pi and phi!=0 and phi!=np.pi
        good_zangle = theta!=np.pi/2

        # intersect with x lines
        righ_dirx_coordx = np.array([signx * self.dx * (0.5 + i)
                                     for i in range(imax[0]) if good_xangle])
        righ_dirx_coordy =  righ_dirx_coordx * np.tan(phi)
        righ_dirx_coordz = RRead.safe_division(self.radius2D(righ_dirx_coordx,
                                                             righ_dirx_coordy),
                                               np.tan(theta))

        left_dirx_coordx = np.array([- signx * self.dx * (0.5 + i)
                                     for i in range(imax[0]) if good_xangle])
        left_dirx_coordy = left_dirx_coordx * np.tan(phi)
        left_dirx_coordz = -RRead.safe_division(self.radius2D(left_dirx_coordx,
                                                              left_dirx_coordy),
                                                np.tan(theta))

        # intersect with y lines
        righ_diry_coordy = np.array([signy * self.dx * (0.5 + i)
                                     for i in range(imax[1]) if good_yangle])
        righ_diry_coordx = RRead.safe_division(righ_diry_coordy, np.tan(phi))
        righ_diry_coordz = RRead.safe_division(self.radius2D(righ_diry_coordx,
                                                             righ_diry_coordy),
                                               np.tan(theta))

        left_diry_coordy = np.array([- signy * self.dx * (0.5 + i)
                                     for i in range(imax[1]) if good_yangle])
        left_diry_coordx = RRead.safe_division(left_diry_coordy,np.tan(phi))
        left_diry_coordz = -RRead.safe_division(self.radius2D(left_diry_coordx,
                                                              left_diry_coordy),
                                                np.tan(theta))

        # intersect with z lines
        righ_dirz_coordz = np.array([signz * self.dx * (0.5 + i)
                                     for i in range(imax[2]) if good_zangle])
        righ_dirz_coordx = righ_dirz_coordz*np.cos(phi)*np.tan(theta)
        righ_dirz_coordy = righ_dirz_coordz*np.sin(phi)*np.tan(theta)

        left_dirz_coordz = np.array([- signz * self.dx * (0.5 + i)
                                     for i in range(imax[2]) if good_zangle])
        left_dirz_coordx = left_dirz_coordz*np.cos(phi)*np.tan(theta)
        left_dirz_coordy = left_dirz_coordz*np.sin(phi)*np.tan(theta)

        #================================================================
        # Calculate distance between points as weighted_dx
        #================================================================
    
        # ===== RIGHT
        # put all the coordinates together
        righ_coordx = [0.0]+list(righ_dirx_coordx)+list(righ_diry_coordx)+list(righ_dirz_coordx)
        righ_coordy = [0.0]+list(righ_dirx_coordy)+list(righ_diry_coordy)+list(righ_dirz_coordy)
        righ_coordz = [0.0]+list(righ_dirx_coordz)+list(righ_diry_coordz)+list(righ_dirz_coordz)
        
        # calc radius of each coordinate
        righ_radius_dirx = [self.radius(xi, yi, zi) 
                            for xi, yi, zi in zip(righ_dirx_coordx, 
                                                  righ_dirx_coordy, 
                                                  righ_dirx_coordz)]
        righ_radius_diry = [self.radius(xi, yi, zi) 
                            for xi, yi, zi in zip(righ_diry_coordx, 
                                                  righ_diry_coordy, 
                                                  righ_diry_coordz)]
        righ_radius_dirz = [self.radius(xi, yi, zi) 
                            for xi, yi, zi in zip(righ_dirz_coordx, 
                                                  righ_dirz_coordy, 
                                                  righ_dirz_coordz)]
        righ_radius = [0.0] + righ_radius_dirx + righ_radius_diry + righ_radius_dirz
        
        # max radius
        righ_max_radius = []
        for radius_dir in [righ_radius_dirx, righ_radius_diry, righ_radius_dirz]:
            if radius_dir!=[]:
                righ_max_radius += [np.max(radius_dir)]
        righ_max_radius = np.min(righ_max_radius)
        if righ_max_radius<rmax:
            righ_max_radius = rmax
        
        # sort the points according to radius size
        righ_coordx = [i for _,i in sorted(zip(righ_radius, righ_coordx))]
        righ_coordy = [i for _,i in sorted(zip(righ_radius, righ_coordy))]
        righ_coordz = [i for _,i in sorted(zip(righ_radius, righ_coordz))]
        righ_radius.sort()
        
        # TO DO: make sure the desired data point is included
        ix = np.argmin(abs(np.array(righ_coordx)-xwanted))
        iy = np.argmin(abs(np.array(righ_coordy)-ywanted))
        iz = np.argmin(abs(np.array(righ_coordz)-zwanted))
        
        # remove repeated radius values and exclude points > max_radius
        righ_coordx_retained, righ_coordy_retained, righ_coordz_retained = [], [], []
        righ_radius_retained = []
        for i, r in enumerate(righ_radius):
            if np.round(r, 3) not in righ_radius_retained and np.round(r, 3)<righ_max_radius:
                righ_radius_retained += [np.round(r, 3)]
                righ_coordx_retained += [righ_coordx[i]]
                righ_coordy_retained += [righ_coordy[i]]
                righ_coordz_retained += [righ_coordz[i]]
        righ_coordx = righ_coordx_retained
        righ_coordy = righ_coordy_retained
        righ_coordz = righ_coordz_retained
        righ_coord  = [righ_coordx, righ_coordy, righ_coordz]
        righ_radius = righ_radius_retained

        # calculate length between each intersecting points to get weight
        righ_weighted_dx = np.array([np.sqrt((righ_coordx[i+1]-righ_coordx[i])**2
                                             + (righ_coordy[i+1]-righ_coordy[i])**2
                                             + (righ_coordz[i+1]-righ_coordz[i])**2) 
                                     for i in range(len(righ_coordx)-1)])
        self.check('1', theta, phi, np.sum(righ_weighted_dx), 
                   self.radius(righ_coordx[-1], righ_coordy[-1], righ_coordz[-1]))

        # ===== LEFT
        # put all the coordinates together
        left_coordx = [0.0]+list(left_dirx_coordx)+list(left_diry_coordx)+list(left_dirz_coordx)
        left_coordy = [0.0]+list(left_dirx_coordy)+list(left_diry_coordy)+list(left_dirz_coordy)
        left_coordz = [0.0]+list(left_dirx_coordz)+list(left_diry_coordz)+list(left_dirz_coordz)
        
        # calc radius of each coordinate
        left_radius_dirx = [self.radius(xi, yi, zi) 
                            for xi, yi, zi in zip(left_dirx_coordx, 
                                                  left_dirx_coordy, 
                                                  left_dirx_coordz)]
        left_radius_diry = [self.radius(xi, yi, zi) 
                            for xi, yi, zi in zip(left_diry_coordx, 
                                                  left_diry_coordy, 
                                                  left_diry_coordz)]
        left_radius_dirz = [self.radius(xi, yi, zi) 
                            for xi, yi, zi in zip(left_dirz_coordx, 
                                                  left_dirz_coordy, 
                                                  left_dirz_coordz)]
        left_radius = [0.0] + left_radius_dirx + left_radius_diry + left_radius_dirz
        
        # max radius
        left_max_radius = []
        for radius_dir in [left_radius_dirx, left_radius_diry, left_radius_dirz]:
            if radius_dir!=[]:
                left_max_radius += [np.max(radius_dir)]
        left_max_radius = np.min(left_max_radius)
        if left_max_radius<rmax:
            left_max_radius = rmax
        
        # sort the points according to radius size
        left_coordx = [i for _,i in sorted(zip(left_radius, left_coordx))]
        left_coordy = [i for _,i in sorted(zip(left_radius, left_coordy))]
        left_coordz = [i for _,i in sorted(zip(left_radius, left_coordz))]
        left_radius.sort()
        
        # remove repeated radius values and exclude points > max_radius
        left_coordx_retained, left_coordy_retained, left_coordz_retained = [], [], []
        left_radius_retained = []
        for i, r in enumerate(left_radius):
            if np.round(r, 3) not in left_radius_retained and np.round(r, 3)<left_max_radius:
                left_radius_retained += [np.round(r, 3)]
                left_coordx_retained += [left_coordx[i]]
                left_coordy_retained += [left_coordy[i]]
                left_coordz_retained += [left_coordz[i]]
        left_coordx = left_coordx_retained
        left_coordy = left_coordy_retained
        left_coordz = left_coordz_retained
        left_coord = [left_coordx, left_coordy, left_coordz]
        left_radius = left_radius_retained
        
        # calculate length between each intersecting points to get weight
        left_weighted_dx = np.array([np.sqrt((left_coordx[i+1]-left_coordx[i])**2
                                             + (left_coordy[i+1]-left_coordy[i])**2
                                             + (left_coordz[i+1]-left_coordz[i])**2) 
                                     for i in range(len(left_coordx)-1)])
        self.check('2', theta, phi, np.sum(left_weighted_dx), 
                   self.radius(left_coordx[-1], left_coordy[-1], left_coordz[-1]))

        #================================================================
        # Get associated index in the grid
        #================================================================

        # ===== RIGHT
        # the average between the line and grid intersections 
        # will be closer to the position of the data points
        # x indice
        righ_av_coordx = [np.average([righ_coordx[i], righ_coordx[i+1]]) 
                          for i in range(len(righ_coordx)-1)]
        righ_ix = [self.closest_idx(i) for i in righ_av_coordx]
        
        # y indice
        righ_av_coordy = [np.average([righ_coordy[i], righ_coordy[i+1]]) 
                          for i in range(len(righ_coordy)-1)]
        righ_iy = [self.closest_idx(i) for i in righ_av_coordy]
        
        # z indice
        righ_av_coordz = [np.average([righ_coordz[i], righ_coordz[i+1]]) 
                          for i in range(len(righ_coordz)-1)]
        righ_iz = [self.closest_idx(i) for i in righ_av_coordz]
        
        # put the indices together and get their corresponding radius
        righ_indices = [(ix, iy, iz) for ix, iy, iz in zip(righ_ix, righ_iy, righ_iz)]
        #if idx_grid_point_considered not in righ_indices:
        #    raise ValueError("The data point I'm after isn't in the indice list")
        righ_radius_from_indices = self.radius_idx(righ_indices, cut=False)
        
        # ===== CHECK if there are gaps between the indices
        gap_between_indices = np.array([abs(np.array(righ_indices[i]) 
                                            - np.array(righ_indices[i+1]))
                                        for i in range(len(righ_indices)-1)])
        # fix because periodic
        for i in range(len(righ_indices)-1):
            for j in range(3):
                if gap_between_indices[i,j] == self.N - 1:
                    gap_between_indices[i,j] = 1
                    
        # take max difference
        gap_between_indices = [np.max(gap_between_indices[i]) 
                               for i in range(len(righ_indices)-1)]
        # if bigger than 1 then error
        give_warning = False
        for i in range(len(righ_indices)-1):
            if gap_between_indices[i] > 1:
                give_warning = True
        if give_warning:
            print()
            print('righ_coordx ', righ_coordx)
            print('righ_av_coordx ', righ_av_coordx)
            print()
            print('righ_coordy ', righ_coordy)
            print('righ_av_coordy ', righ_av_coordy)
            print()
            print('righ_coordz ', righ_coordz)
            print('righ_av_coordz ', righ_av_coordz)
            print()
            print('righ_indices ', righ_indices)
            print('gap_between_indices ', gap_between_indices)
            raise ValueError('Right, gap between indices too big')
        
        # ===== LEFT
        # x indice
        left_av_coordx = [np.average([left_coordx[i], left_coordx[i+1]]) 
                          for i in range(len(left_coordx)-1)]
        left_ix = [self.closest_idx(i) for i in left_av_coordx]
        
        # y indice
        left_av_coordy = [np.average([left_coordy[i], left_coordy[i+1]]) 
                          for i in range(len(left_coordy)-1)]
        left_iy = [self.closest_idx(i) or i in left_av_coordy]
        
        # z indice
        left_av_coordz = [np.average([left_coordz[i], left_coordz[i+1]]) 
                          for i in range(len(left_coordz)-1)]
        left_iz = [self.closest_idx(i) for i in left_av_coordz]
        
        # put the indices together and get their corresponding radius
        left_indices = [(ix, iy, iz) for ix, iy, iz in zip(left_ix, left_iy, left_iz)]
        #if idx_grid_point_considered not in left_indices:
        #    raise ValueError("Left, the data point I'm after isn't in the indice list")
        left_radius_from_indices = self.radius_idx(left_indices)
        
        # CHECK if there are gaps between the indices
        gap_between_indices = np.array([abs(np.array(left_indices[i]) 
                                            - np.array(left_indices[i+1]))
                                        for i in range(len(left_indices)-1)])
        for i in range(len(left_indices)-1):
            for j in range(3):
                if gap_between_indices[i,j] == self.N - 1:
                    gap_between_indices[i,j] = 1
        gap_between_indices = [np.max(gap_between_indices[i]) 
                               for i in range(len(left_indices)-1)]
        give_warning = False
        for i in range(len(left_indices)-1):
            if gap_between_indices[i] > 1:
                give_warning = True
        if give_warning:
            raise ValueError('Left, gap between indices too big')
        
        #================================================================
        # Get associated index in the grid
        #================================================================
        left_values = left_coord, left_radius, left_indices, left_radius_from_indices, left_weighted_dx
        righ_values = righ_coord, righ_radius, righ_indices, righ_radius_from_indices, righ_weighted_dx
        return left_values, righ_values

    def lin_fit_zero(self, f, r):
        if r[0]==r[1]:
            return np.array(r[0])
        else:
            a = RRead.safe_division(f[0]-f[1], r[0]-r[1])
            b = f[0]-r[0]*a
            return RRead.safe_division(-b, a)
        
    def lin_fit(self, f, r, x):
        a = RRead.safe_division((f[0]-f[1]), (r[0]-r[1]))
        b = f[0]-r[0]*a
        return a*x+b
    
    def tri_lin_fit(self, K, idx0, bc0, idx1, X):
        x,y,z = X
        x0, y0, z0 = list(idx0)
        x0_phy = self.xyz_periodic(x0, boundary_crossed=bc0[0])[0]
        y0_phy = self.xyz_periodic(y0, boundary_crossed=bc0[1])[0]
        z0_phy = self.xyz_periodic(z0, boundary_crossed=bc0[2])[0]
        x1, y1, z1 = list(idx1)
        
        if x0==x1:
            if x0_phy<x:
                if x1==self.N-1:
                    x1=0
                else:
                    x1+=1
            else:
                x1+= -1
        if y0==y1:
            if y0_phy<y:
                if y1==self.N-1:
                    y1=0
                else:
                    y1+=1
            else:
                y1+= -1
        if z0==z1:
            if z0_phy<z:
                if z1==self.N-1:
                    z1=0
                else:
                    z1+=1
            else:
                z1+= -1
            
        x1_phy = self.xyz_periodic(x1, boundary_crossed=(bc0[0] or self.xyz_periodic(x1)[1]))[0]
        y1_phy = self.xyz_periodic(y1, boundary_crossed=(bc0[1] or self.xyz_periodic(y1)[1]))[0]
        z1_phy = self.xyz_periodic(z1, boundary_crossed=(bc0[2] or self.xyz_periodic(z1)[1]))[0]
        
        K_y0z0 = self.lin_fit([K[x0,y0,z0], K[x1,y0,z0]], [x0_phy,x1_phy], x)
        K_y1z0 = self.lin_fit([K[x0,y1,z0], K[x1,y1,z0]], [x0_phy,x1_phy], x)
        K_z0 = self.lin_fit([K_y0z0, K_y1z0], [y0_phy,y1_phy], y)
        
        K_y0z1 = self.lin_fit([K[x0,y0,z1], K[x1,y0,z1]], [x0_phy,x1_phy], x)
        K_y1z1 = self.lin_fit([K[x0,y1,z1], K[x1,y1,z1]], [x0_phy,x1_phy], x)
        K_z1 = self.lin_fit([K_y0z1, K_y1z1], [y0_phy,y1_phy], y)
        return self.lin_fit([K_z0, K_z1], [z0_phy,z1_phy], z)
    
    def extrapolate(self, idx_location, indices, radius_from_indices, 
                    coord, radius, f_in_grid, smaller_than_r=True):
        idx_location_in_grid = 0
        while radius_from_indices[idx_location_in_grid]<radius[idx_location]: idx_location_in_grid+=1
        if smaller_than_r:
            idx_location_in_grid-=1
        #idx_location_in_grid = np.argmin(abs(radius_from_indices-radius[idx_location]))
        
        coordx, coordy, coordz = coord
        boundary_crossed_x = self.xyz_periodic([i[0] for i in indices[:idx_location_in_grid+1]])[1]
        boundary_crossed_y = self.xyz_periodic([i[1] for i in indices[:idx_location_in_grid+1]])[1]
        boundary_crossed_z = self.xyz_periodic([i[2] for i in indices[:idx_location_in_grid+1]])[1]
        bc0 = [boundary_crossed_x, boundary_crossed_y, boundary_crossed_z]
        i = 1
        while indices[idx_location_in_grid]==indices[idx_location_in_grid+i]: i+=1
        
        f_extrap = self.tri_lin_fit(f_in_grid, indices[idx_location_in_grid], 
                                    bc0, indices[idx_location_in_grid+i], 
                                    [coordx[idx_location], coordy[idx_location], coordz[idx_location]])
        return f_extrap
    
    def get_TA_Radius(self, gdet_in_grid, K_in_grid, theta, phi, rmaxfac):
        Radius_comoving = []
        Radius_proper = []
        for ti, pi in zip(theta, phi): 
            indices, weighted_dx = self.get_weighted_indices_forTA(ti, pi, rmaxfac, K_in_grid)
            for i, w in zip(indices, weighted_dx):
                gdet  = np.array([gdet_in_grid[idx] for idx in i])
                Radius_comoving += [np.sum(w)]
                Radius_proper += [np.sum(w*gdet**(1/6))]
            
        return np.average(Radius_comoving),  np.average(Radius_proper)
    
    def get_region_Radius(self, radius, gdet_in_grid, region, theta, phi, rmaxfac):
        Radius_comoving = []
        Radius_proper = []
        for ti, pi in zip(theta, phi): 
            indices, weighted_dx = self.get_weighted_indices_for_region(radius, region, ti, pi, rmaxfac)
            for i, w in zip(indices, weighted_dx):
                gdet  = np.array([gdet_in_grid[idx] for idx in i])
                Radius_comoving += [np.sum(w)]
                Radius_proper += [np.sum(w*gdet**(1/6))]
        return np.average(Radius_comoving),  np.average(Radius_proper)
    
    def get_proper_Radius_of_grid_point(self, gdet_in_grid, idx_grid_point_considered):
        if idx_grid_point_considered == (self.ixOD, self.iyOD, self.izOD):
            return 0.0, 0.0, 0.0
        else:
            values = self.get_weighted_indices_for_gridpoint(idx_grid_point_considered)
            indices, weighted_dx, theta, phi = values

            gdet  = np.array([gdet_in_grid[idx] for idx in indices])
            Radius_proper = np.sum(weighted_dx*gdet**(1/6))
            return Radius_proper, theta, phi
    
    def get_Radius_of_grid_point(self, idx_grid_point_considered):
        if idx_grid_point_considered == (self.ixOD, self.iyOD, self.izOD):
            return 0.0, 0.0, 0.0
        else:
            values = self.get_weighted_indices_for_gridpoint(idx_grid_point_considered)
            indices, weighted_dx, theta, phi = values
            Radius = np.sum(weighted_dx)
            return Radius, theta, phi
    
    def get_weighted_indices_forTA(self, theta, phi, rmaxfac, K_in_grid):
        left_values, righ_values = self.idx_along_r(theta, phi, rmaxfac)
        final_indices, final_weighted_dx = [], []
        for values in [left_values, righ_values]:
            coord, radius, indices, radius_from_indices, weighted_dx = values
            #========================================
            # Extrapolate the last point
            #========================================

            # ======= 1st Get the K at the intermediary grid points around TA
            # == Determine intermediary grid points around TA
            K     = np.array([K_in_grid[idx] for idx in indices])
            pre_TA = 0
            while np.sign(K[0])==np.sign(K[pre_TA]): pre_TA +=1
            pre_TA -=1
            self.check('3', theta, phi, np.sum(weighted_dx[:pre_TA]), 
                       self.radiusOD(coord[0][pre_TA], coord[1][pre_TA], coord[2][pre_TA]))


            # == Trilinear extrapolation to get K
            K_extrap_pre_TA = self.extrapolate(pre_TA, indices, radius_from_indices, 
                                               coord, radius, K_in_grid)
            while K_extrap_pre_TA<0:
                pre_TA-=1
                K_extrap_pre_TA = self.extrapolate(pre_TA, indices, radius_from_indices, 
                                                   coord, radius, K_in_grid)
            K_extrap_prenext_TA = self.extrapolate(pre_TA+1, indices, radius_from_indices, 
                                                   coord, radius, K_in_grid)
            while K_extrap_prenext_TA>0:
                pre_TA+=1
                K_extrap_pre_TA = self.extrapolate(pre_TA, indices, radius_from_indices, 
                                                   coord, radius, K_in_grid)
                K_extrap_prenext_TA = self.extrapolate(pre_TA+1, indices, radius_from_indices, 
                                                       coord, radius, K_in_grid)


            post_TA = pre_TA+1
            while np.round(radius[pre_TA],1)==np.round(radius[post_TA],1): post_TA+=1
            K_extrap_post_TA = self.extrapolate(post_TA, indices, radius_from_indices, 
                                                coord, radius, K_in_grid, 
                                                smaller_than_r=False)
            while K_extrap_post_TA>0:
                post_TA+=1
                K_extrap_post_TA = self.extrapolate(post_TA, indices, radius_from_indices, 
                                                    coord, radius, K_in_grid, 
                                                    smaller_than_r=False)

            if np.sign(K_extrap_pre_TA)==np.sign(K_extrap_post_TA):
                print('WARNING: sign pre and post extrapolated K are the same, both:', 
                      np.sign(K_extrap_pre_TA))

            # == Extrapolation to TA location
            TA_coordx = self.lin_fit_zero([K_extrap_pre_TA, K_extrap_post_TA], 
                                          [coord[0][pre_TA], coord[0][post_TA]])
            TA_coordy = self.lin_fit_zero([K_extrap_pre_TA, K_extrap_post_TA], 
                                          [coord[1][pre_TA], coord[1][post_TA]])
            TA_coordz = self.lin_fit_zero([K_extrap_pre_TA, K_extrap_post_TA], 
                                          [coord[2][pre_TA], coord[2][post_TA]])
            radius_TA = self.lin_fit_zero([K_extrap_pre_TA, K_extrap_post_TA], 
                                          [radius[pre_TA], radius[post_TA]])

            # == Get indice 
            TA_indices = (np.argmin(abs(self.Lin.d1x-TA_coordx)), 
                          np.argmin(abs(self.Lin.d1x-TA_coordy)), 
                          np.argmin(abs(self.Lin.d1x-TA_coordz)))
            final_indices += [[indices[i] for i in range(pre_TA)]+[TA_indices]]

            # == Get weight        
            last_weighted_dx = np.sqrt((coord[0][pre_TA]-TA_coordx)**2
                                       +(coord[1][pre_TA]-TA_coordy)**2
                                       +(coord[2][pre_TA]-TA_coordz)**2)
            final_weighted_dx += [np.append(weighted_dx[:pre_TA], last_weighted_dx)]
            self.check('6', theta, phi, np.sum(final_weighted_dx[-1]), 
                       self.radiusOD(TA_coordx, TA_coordy, TA_coordz))
        
        return final_indices, final_weighted_dx
    
    
    def get_weighted_indices_for_region(self, wanted_radius, region, theta, phi, rmaxfac):
        left_values, righ_values = self.idx_along_r(theta, phi, rmaxfac)
        final_indices, final_weighted_dx = [], []
        for values in [left_values, righ_values]:
            coord, radius, indices, radius_from_indices, weighted_dx = values
            
            i = 0
            ind_in_region = []
            while indices[i] in region:
                ind_in_region += [indices]
                i+=1
            
            pre_final_indices = indices[:i]
            pre_final_weighted_dx = list(weighted_dx[:i])
            if wanted_radius<np.sum(pre_final_weighted_dx):
                pre_final_weighted_dx += [wanted_radius-radius[i]]
                pre_final_indices += [indices[i-1]]
            elif wanted_radius>np.sum(pre_final_weighted_dx):
                pre_final_weighted_dx += [wanted_radius-radius[i]]
                pre_final_indices += [indices[i]]
                
            if abs(wanted_radius/np.sum(pre_final_weighted_dx)-1)>1e-3:
                print('WARNING:', wanted_radius, np.sum(pre_final_weighted_dx))
            final_indices += [pre_final_indices]
            final_weighted_dx += [pre_final_weighted_dx]
        
        return final_indices, final_weighted_dx
    
    def get_weighted_indices_for_gridpoint(self, idx_grid_point_considered):
        # Cartesian coordinates
        # Shift around OD
        x = -self.xOD + self.Lin.d1x[idx_grid_point_considered[0]]
        y = -self.yOD + self.Lin.d1x[idx_grid_point_considered[1]]
        z = -self.zOD + self.Lin.d1x[idx_grid_point_considered[2]]
        if x>=self.L/2:
            x = x - self.L
        if y>=self.L/2:
            y = y - self.L
        if z>=self.L/2:
            z = z - self.L

        # Spherical coordinates
        # radius
        r = self.radius(x, y, z)
        # inclination
        theta = np.arccos(RRead.safe_division(z, r))
        if r==0:
            theta = 0
        # azimuth
        phi = np.arccos(RRead.safe_division(x, self.radius2D(x, y)))
        if y<0:
            phi = 2*np.pi - phi
        if x==0 and y==0:
            phi = 0
            
        xc = r * np.cos(phi) * np.sin(theta)
        yc = r * np.sin(phi) * np.sin(theta)
        zc = r * np.cos(theta)
        self.check('x', theta, phi, x, xc)
        self.check('y', theta, phi, y, yc)
        self.check('z', theta, phi, z, zc)
            
        buffer = 4
        rmax = r + np.sqrt(3) * self.dx
        nr = r + buffer * np.sqrt(3) * self.dx
        nnr = nr + buffer * np.sqrt(3) * self.dx
        rx, ry, rz = [i for _,i in sorted(zip(np.array([x,y,z]), 
                                              np.array([r,nr,nnr])))]
        xmax = rx * np.cos(phi) * np.sin(theta)
        ymax = ry * np.sin(phi) * np.sin(theta)
        zmax = ry * np.cos(theta)
        imax = [int(abs(i/self.dx))+buffer for i in [xmax, ymax, zmax]]
            
        left_values, righ_values = self.idx_along_r(theta, phi, imax, rmax, 
                                                    idx_grid_point_considered, x, y, z)
        coord, radius, indices, radius_from_indices, weighted_dx = righ_values
        coordL, radiusL, indicesL, radius_from_indicesL, weighted_dxL = left_values
        
        # shift everything back around OD
        fcoord = [coord[0]+self.xOD, coord[1]+self.yOD, coord[2]+self.zOD]
        findices = []
        for idx in indices:
            ix, iy, iz = idx
            ix = ix - self.ixOD
            iy = iy - self.iyOD
            iz = iz - self.izOD
            if ix<0:
                ix = self.N + ix
            if iy<0:
                iy = self.N + iy
            if iz<0:
                iz = self.N + iz
            if ix>self.N-1:
                ix = ix - self.N
            if iy>self.N-1:
                iy = iy - self.N
            if iz>self.N-1:
                iz = iz - self.N
            findices += [(ix, iy, iz)]
            
        # cut data to wanted data point
        try:
            i = findices.index(idx_grid_point_considered)+1
        except:
            print(x, y, z)
            print(indices)
            print(indicesL)
            print(findices)
            i = findices.index(idx_grid_point_considered)+1
        
        final_indices = findices[:i]
        final_weighted_dx = list(weighted_dx[:i])
        if r<np.sum(final_weighted_dx):
            final_weighted_dx += [r-radius[i]]
            final_indices += [findices[i-1]]
        elif r>np.sum(final_weighted_dx):
            final_weighted_dx += [r-radius[i]]
            final_indices += [findices[i]]
        
        return final_indices, final_weighted_dx, theta, phi
    
    def get_delta_in_loc(self, theta, phi, rmaxfac, delta_in_grid, wanted_length):
        left_values, righ_values = self.idx_along_r(theta, phi, rmaxfac)
        
        delta_saved = []
        for values in [left_values, righ_values]:
            coord, radius, indices, radius_from_indices, weighted_dx = values
            pre_loc = 0
            while radius[pre_loc]<wanted_length: pre_loc+=1
            pre_loc-=1
            
            delta_extrap_pre_loc = self.extrapolate(pre_loc, indices, radius_from_indices, 
                                                    coord, radius, delta_in_grid)
            delta_extrap_post_loc = self.extrapolate(pre_loc+1, indices, radius_from_indices, 
                                                     coord, radius, delta_in_grid, 
                                                     smaller_than_r=False)
            delta_loc = self.lin_fit([delta_extrap_pre_loc, delta_extrap_post_loc], 
                                     [radius[pre_loc], radius[pre_loc+1]], wanted_length)
            delta_saved += [delta_loc]
        return delta_saved
    
    def check(self, i, theta, phi, A, B):
        if B!=0 and A!=0:
            if abs(A/B - 1)>1e-8:
                print('WARNING:' + i 
                      + ' theta={:.2f}, phi={:.2f}, {:.5f}!={:.5f}'.format(theta*180/np.pi, 
                                                                           phi*180/np.pi, A, B), ' ')

            
            
            
            
            
            
            
            
        