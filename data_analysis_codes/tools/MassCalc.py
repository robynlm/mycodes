import numpy as np

class MassCalcClass:
    def __init__(self, param, Lin, OD_Loc):
        self.dx = param['dx']
        self.L = param['Lx']
        self.N = param['Nx']
        self.Lin = Lin
        self.xOD = self.Lin.d3xyz[OD_Loc[0][0]]
        self.yOD = self.Lin.d3xyz[OD_Loc[1][0]]
        self.zOD = self.Lin.d3xyz[OD_Loc[2][0]]
        self.xyzOD = np.array([self.xOD, self.yOD, self.zOD])
        
    def CalcMass(self, rho, gdet, r):
        weight = self.getWeightsOfWholeBox(r)
        return np.sum(rho*weight*np.sqrt(gdet))*(self.dx**3), weight
    
    def Average_in_Sphere(self, phi, gdet, r):
        weights = self.getWeightsOfWholeBox(r)
        return np.sum(phi*weights*np.sqrt(gdet))/np.sum(weights*np.sqrt(gdet))
    
    def getvar(self, rho, gdet, r):
        weights = self.getWeightsOfWholeBox(r)
        delta = rho/self.Lin.evo.rho(self.Lin.t_initial) - 1
        Volume = weights*np.sqrt(gdet)
        Mass = np.sum(rho*Volume)*(self.dx**3)
        delta_av = np.sum(delta*Volume)/np.sum(Volume)
        return [r, Mass, Mass*self.Lin.Massfac, delta_av]
    
    def Collect_points(self,r):
        P = []
        for ix in range(self.N):
            for iy in range(self.N):
                for iz in range(self.N):
                    pos1, pos2, pos3, pos4, pos5, pos6, pos7, pos8 = self.get8positions(self.Lin.d3xyz[ix], self.Lin.d3xyz[iy], self.Lin.d3xyz[iz])
                    pos1check = self.check_if_in_radius(pos1, r)
                    pos2check = self.check_if_in_radius(pos2, r)
                    pos3check = self.check_if_in_radius(pos3, r)
                    pos4check = self.check_if_in_radius(pos4, r)
                    pos5check = self.check_if_in_radius(pos5, r)
                    pos6check = self.check_if_in_radius(pos6, r)
                    pos7check = self.check_if_in_radius(pos7, r)
                    pos8check = self.check_if_in_radius(pos8, r)
                    poscheckall = np.array([pos1check, pos2check, pos3check, pos4check, pos5check, pos6check, pos7check, pos8check])
                    poscheck = np.sum(poscheckall)
                    
                    if poscheck>0:
                        nexttopos1 = [pos2, pos3, pos7]
                        nexttopos2 = [pos1, pos4, pos8]
                        nexttopos3 = [pos1, pos4, pos5]
                        nexttopos4 = [pos2, pos3, pos6]
                        nexttopos5 = [pos3, pos6, pos7]
                        nexttopos6 = [pos4, pos5, pos8]
                        nexttopos7 = [pos1, pos5, pos8]
                        nexttopos8 = [pos2, pos6, pos7]
                        posall = np.array([pos1, pos2, pos3, pos4, pos5, pos6, pos7, pos8])
                        nexttoposall = np.array([nexttopos1, nexttopos2, nexttopos3, nexttopos4, nexttopos5, nexttopos6, nexttopos7, nexttopos8])
                        if poscheck==7:
                            P += self.Points1point(r, posall[np.where(~poscheckall)], nexttoposall[np.where(~poscheckall)])
                        elif poscheck==6:
                            P += self.Points2points(r, posall[np.where(~poscheckall)], nexttoposall[np.where(~poscheckall)])
                        elif poscheck==5:
                            P += self.Points3points(r, posall[np.where(~poscheckall)], nexttoposall[np.where(~poscheckall)])
                        elif poscheck==4:
                            P += self.Points4points(r, posall[np.where(poscheckall)], nexttoposall[np.where(poscheckall)])
                        elif poscheck==3:
                            P += self.Points3points(r, posall[np.where(poscheckall)], nexttoposall[np.where(poscheckall)])
                        elif poscheck==2:
                            P += self.Points2points(r, posall[np.where(poscheckall)], nexttoposall[np.where(poscheckall)])
                        elif poscheck==1:
                            P += self.Points1point(r, posall[np.where(poscheckall)], nexttoposall[np.where(poscheckall)])
        return P
    
    def getWeightsOfWholeBox(self, r):
        if r>= self.L/2:
            print('ERROR: that radius is too big for me sorry')
        else:
            weight = np.zeros((self.N, self.N, self.N))
            for ix in range(self.N):
                for iy in range(self.N):
                    for iz in range(self.N):
                        weight[ix,iy,iz] = self.getWeights(self.Lin.d3xyz[ix], self.Lin.d3xyz[iy], self.Lin.d3xyz[iz], r)
            return weight
    
    def getWeights(self,x,y,z, r):
        if ~self.check_if_in_radius(np.array([x, y, z]), r-self.dx*np.sqrt(3)) and self.check_if_in_radius(np.array([x, y, z]), r+self.dx*np.sqrt(3)):       
            pos1, pos2, pos3, pos4, pos5, pos6, pos7, pos8 = self.get8positions(x, y, z)
            pos1check = self.check_if_in_radius(pos1, r)
            pos2check = self.check_if_in_radius(pos2, r)
            pos3check = self.check_if_in_radius(pos3, r)
            pos4check = self.check_if_in_radius(pos4, r)
            pos5check = self.check_if_in_radius(pos5, r)
            pos6check = self.check_if_in_radius(pos6, r)
            pos7check = self.check_if_in_radius(pos7, r)
            pos8check = self.check_if_in_radius(pos8, r)
            poscheckall = np.array([pos1check, pos2check, pos3check, pos4check, pos5check, pos6check, pos7check, pos8check])
            poscheck = np.sum(poscheckall)
            
            if poscheck>0:
                nexttopos1 = [pos2, pos3, pos7]
                nexttopos2 = [pos1, pos4, pos8]
                nexttopos3 = [pos1, pos4, pos5]
                nexttopos4 = [pos2, pos3, pos6]
                nexttopos5 = [pos3, pos6, pos7]
                nexttopos6 = [pos4, pos5, pos8]
                nexttopos7 = [pos1, pos5, pos8]
                nexttopos8 = [pos2, pos6, pos7]

                posall = np.array([pos1, pos2, pos3, pos4, pos5, pos6, pos7, pos8])
                nexttoposall = np.array([nexttopos1, nexttopos2, nexttopos3, nexttopos4, nexttopos5, nexttopos6, nexttopos7, nexttopos8])

                Vbox = self.dx**3
                if poscheck==8:
                    return 1.0
                elif poscheck==7:
                    V = self.Volume1point(r, posall[np.where(~poscheckall)], nexttoposall[np.where(~poscheckall)])
                    return 1-V/Vbox
                elif poscheck==6:
                    V = self.Volume2points(r, posall[np.where(~poscheckall)], nexttoposall[np.where(~poscheckall)])
                    return 1-V/Vbox
                elif poscheck==5:
                    V = self.Volume3points(r, posall[np.where(~poscheckall)], nexttoposall[np.where(~poscheckall)])
                    return 1-V/Vbox
                elif poscheck==4:
                    V = self.Volume4points(r, posall[np.where(poscheckall)], nexttoposall[np.where(poscheckall)])
                    return V/Vbox
                elif poscheck==3:
                    V = self.Volume3points(r, posall[np.where(poscheckall)], nexttoposall[np.where(poscheckall)])
                    return V/Vbox
                elif poscheck==2:
                    V = self.Volume2points(r, posall[np.where(poscheckall)], nexttoposall[np.where(poscheckall)])
                    return V/Vbox
                elif poscheck==1:
                    V = self.Volume1point(r, posall[np.where(poscheckall)], nexttoposall[np.where(poscheckall)])
                    return V/Vbox
            else:
                return 0.0
        elif self.check_if_in_radius(np.array([x, y, z]), r-self.dx*np.sqrt(3)):
            return 1.0
        else:
            return 0.0
        
    def get8positions(self,x,y,z):
        dxmax = self.dx/2
        pos1 = np.array([x+dxmax, y+dxmax, z+dxmax])
        pos2 = np.array([x-dxmax, y+dxmax, z+dxmax])
        pos3 = np.array([x+dxmax, y-dxmax, z+dxmax])
        pos4 = np.array([x-dxmax, y-dxmax, z+dxmax])
        pos5 = np.array([x+dxmax, y-dxmax, z-dxmax])
        pos6 = np.array([x-dxmax, y-dxmax, z-dxmax])
        pos7 = np.array([x+dxmax, y+dxmax, z-dxmax])
        pos8 = np.array([x-dxmax, y+dxmax, z-dxmax])
        return pos1, pos2, pos3, pos4, pos5, pos6, pos7, pos8          
    
    def check_if_in_radius(self,pos,r):
        if r<=self.L/4:
            return np.sqrt(np.sum((pos-self.xyzOD)**2)) <= r
        else: 
            if np.sqrt(np.sum((pos-self.xyzOD)**2))<=self.L/4:
                return np.sqrt(np.sum((pos-self.xyzOD)**2)) <= r
            else:# need to use periodic boundaries
                newpos = []
                for coord in range(3):
                    if pos[coord]>self.L/4:
                        newpos += [pos[coord] - self.L]
                    else:
                        newpos += [pos[coord]]
                return np.sqrt(np.sum((np.array(newpos)-self.xyzOD)**2)) <= r
    
    def interpolate_radius_get_distance(self, r, point, neighbourpoint):
        distances_to_point = []
        ixyz = []
            
        periodic_boundaries_are_needed = np.sqrt(np.sum((point-self.xyzOD)**2))>self.L/4
        for npoint in neighbourpoint:
            periodic_boundaries_are_needed = periodic_boundaries_are_needed and np.sqrt(np.sum((npoint-self.xyzOD)**2))>self.L/4
            
        #periodic_boundaries_are_partialy_needed = np.sqrt(np.sum((point-self.xyzOD)**2))<=self.Lin.p.L/4
        #for npoint in neighbourpoint:
        #    periodic_boundaries_are_partialy_needed = periodic_boundaries_are_partialy_needed or np.sqrt(np.sum((npoint-self.xyzOD)**2))<=self.Lin.p.L/4
            
        if periodic_boundaries_are_needed:
            for coord in range(3):
                if point[coord]>self.L/4:
                    point[coord] -= self.L
                for p in neighbourpoint:
                    if p[coord]>self.L/4:
                        p[coord] -= self.L
            
        for p in neighbourpoint:
            interpolated_point = p.copy() #point on the surface between the point in the sphere and it's neighbour
            coordtochange = np.where(point-p!=0)
            direction = np.sign(1 if point[coordtochange]>self.xyzOD[coordtochange] else -1)
            interpolated_point[coordtochange] = direction * np.sqrt( r**2 - np.sum((point-self.xyzOD)**2) + (point[coordtochange]-self.xyzOD[coordtochange])**2 ) + self.xyzOD[coordtochange]
            d_to_point = np.sqrt(np.sum((point-interpolated_point)**2))
            distances_to_point += [d_to_point]
            ixyz += [coordtochange]
            if d_to_point>self.dx:
                print('WARNING: depth > dx')
        return np.array(distances_to_point), ixyz
            
    
    def interpolate_radius(self, r, point, neighbourpoint):
        points = []
        for p in neighbourpoint:
            interpolated_point = p.copy()
            coordtochange = np.where(point-p!=0)
            direction = np.sign(1 if point[coordtochange]>self.xyzOD[coordtochange] else -1)
            interpolated_point[coordtochange] = direction * np.sqrt( r**2 - np.sum((point-self.xyzOD)**2) + (point[coordtochange]-self.xyzOD[coordtochange])**2 ) + self.xyzOD[coordtochange]
            points += [interpolated_point]
        return points
    
    def Volume1point(self, r, point, pointneighbour):  # Half Prism, with a right triangle base
        PrismSides, ixyz = self.interpolate_radius_get_distance(r, point[0], pointneighbour[0])
        # The sides correspond to those connected to the right angle
        V = np.prod(PrismSides) / 6  # 1/2 bcs prism and 1/3 bcs we want pyramid
        if V> (1/6) * self.dx**3:
            print('WARNING: Volume 1 point too big')
        return V
        
        
    def Volume2points(self, r, points, pointneighbours):  # Prism, with a right triangle base
        point1 = points[0]
        point2 = points[1]
        pointneighbour1 = [neighbour for neighbour in pointneighbours[0] if list(neighbour)!=list(point2)]
        pointneighbour2 = [neighbour for neighbour in pointneighbours[1] if list(neighbour)!=list(point1)]
        
        # Prism with point1 at right triangle base
        TriangleBaseSides1, ixyz1 = self.interpolate_radius_get_distance(r, point1, pointneighbour1)
        V1 = self.dx * np.prod(TriangleBaseSides1) / 2
        
        # Prism with point2 at right triangle base
        TriangleBaseSides2, ixyz2 = self.interpolate_radius_get_distance(r, point2, pointneighbour2)
        V2 = self.dx * np.prod(TriangleBaseSides2) / 2
        
        if abs(V1/V2-1)>1e-2: # base triangles are different
            # To find this: Volume = Big half prism - Little half prism
            #                    our extended shape - the extended part
            # I'm making sure the depth/width of the triangles overlap
            dict1 = {ixyz1[i][0][0]:TriangleBaseSides1[i] for i in range(len(ixyz1))}
            dict2 = {ixyz2[i][0][0]:TriangleBaseSides2[i] for i in range(len(ixyz2))}
            idict = list(dict1.keys())
            if V1>V2: #point1 has the bigger base
                BigTriangleWidth = dict1[idict[0]]
                BigTriangleDepth = dict1[idict[1]]
                LitTriangleWidth = dict2[idict[0]]
                LitTriangleDepth = dict2[idict[1]]
            else:     #point2 has the bigger base
                BigTriangleWidth = dict2[idict[0]]
                BigTriangleDepth = dict2[idict[1]]
                LitTriangleWidth = dict1[idict[0]]
                LitTriangleDepth = dict1[idict[1]]
                
            LitPrismHeight = LitTriangleWidth * self.dx / ( BigTriangleWidth - LitTriangleWidth )
            VBigPrism = BigTriangleWidth * BigTriangleDepth * (LitPrismHeight + self.dx)/ 2
            VLitPrism = LitTriangleWidth * LitTriangleDepth * LitPrismHeight / 2
            V = (VBigPrism - VLitPrism) / 3
        else:
            V = np.average([V1, V2])
            
        if V>(1/2)*self.dx**3:
            print('WARNING: Volume 2 points too big')
        return V
    
    def Volume3points(self, r, points, pointneighbours):
        # This is like the heptahedron case for the 4 point case
        # But c1 = c2 and the box width is smaller = dist(A, C)
        point1 = points[0]
        point2 = points[1]
        point3 = points[2]
        pointneighbour1 = [neighbour for neighbour in pointneighbours[0] if list(neighbour)!=list(point2) and list(neighbour)!=list(point3)]
        pointneighbour2 = [neighbour for neighbour in pointneighbours[1] if list(neighbour)!=list(point1) and list(neighbour)!=list(point3)]
        pointneighbour3 = [neighbour for neighbour in pointneighbours[2] if list(neighbour)!=list(point1) and list(neighbour)!=list(point2)]
        
        depth1, ixyz1 = self.interpolate_radius_get_distance(r, point1, pointneighbour1)
        depth2, ixyz2 = self.interpolate_radius_get_distance(r, point2, pointneighbour2)
        depth3, ixyz3 = self.interpolate_radius_get_distance(r, point3, pointneighbour3)
        depth = [depth1, depth2, depth3]
        ixyz = [ixyz1, ixyz2, ixyz3]
        
        # ID points
        for i, d in enumerate(depth):
            if len(d)==1:
                iA = i
        iB, iD = np.delete(np.arange(3), iA) #index of other points
        ixyzA = ixyz[iA]
        ixyzB = ixyz[iB]
        ixyzD = ixyz[iD]
        depthA = depth[iA]
        depthB = depth[iB]
        depthD = depth[iD]
            
        # ID depths
        a1 = depthA
        dictB = {ixyzB[i][0][0]:depthB[i] for i in range(len(ixyzB))}
        dictD = {ixyzD[i][0][0]:depthD[i] for i in range(len(ixyzD))}
        for keyB in dictB:
            for keyD in dictD:
                if keyB==keyD:
                    b1 = dictB[keyB]
                    d1 = dictD[keyD]
                    takenkey = keyB
        for keyB in dictB:
            if keyB!=takenkey:
                b2 = dictB[keyB]
        for keyD in dictD:
            if keyD!=takenkey:
                d2 = dictD[keyD]
        
        Bleg = b1 * self.dx / (a1 - b1)
        Dleg = d1 * self.dx / (a1 - d1)
            
        BigPyramid = (self.dx+Bleg) * (a1) * (self.dx+Dleg) / 6
        BPyramid = Bleg * b1 * b2 / 6
        DPyramid = Dleg * d1 * d2 / 6
            
        V = BigPyramid - BPyramid - DPyramid
        if V>(5/6) * self.dx**3:
                print('WARNING: Volume 3 points too big')
        return V
    
    def Volume4points(self, r, points, pointneighbours):
        point1 = points[0]
        point2 = points[1]
        point3 = points[2]
        point4 = points[3]
        
        pointneighbour1 = [neighbour for neighbour in pointneighbours[0] if list(neighbour)!=list(point2) and list(neighbour)!=list(point3) and list(neighbour)!=list(point4)]
        pointneighbour2 = [neighbour for neighbour in pointneighbours[1] if list(neighbour)!=list(point1) and list(neighbour)!=list(point3) and list(neighbour)!=list(point4)]
        pointneighbour3 = [neighbour for neighbour in pointneighbours[2] if list(neighbour)!=list(point1) and list(neighbour)!=list(point2) and list(neighbour)!=list(point4)]
        pointneighbour4 = [neighbour for neighbour in pointneighbours[3] if list(neighbour)!=list(point1) and list(neighbour)!=list(point2) and list(neighbour)!=list(point3)]
        
        depth1, ixyz1 = self.interpolate_radius_get_distance(r, point1, pointneighbour1)
        depth2, ixyz2 = self.interpolate_radius_get_distance(r, point2, pointneighbour2)
        depth3, ixyz3 = self.interpolate_radius_get_distance(r, point3, pointneighbour3)
        depth4, ixyz4 = self.interpolate_radius_get_distance(r, point4, pointneighbour4)
        depth = [depth1, depth2, depth3, depth4]
        if len(depth1)+len(depth2)+len(depth3)+len(depth4)==4:
            # irregular prism: Volume = Volume Pave (height min(depth)) 
            #                         + Volume Half Pave (height max(depth)-min(depth))
            V = ( np.max(depth) + np.min(depth) ) * self.dx * self.dx / 2
            if V > self.dx**3:
                print('WARNING: Volume 4 points, 1st case, too big')
        else: # Outch, this is a heptahedron
            # Which one is A? the one where len(depth)==0
            for i, d in enumerate(depth):
                if len(d)==0:
                    iA = i
            # Which one is B, C, D ?
            ixyz = [ixyz1, ixyz2, ixyz3, ixyz4]
            iB, iC, iD = np.delete(np.arange(4), iA) #index of other points
            ixyzB = ixyz[iB]
            ixyzC = ixyz[iC]
            ixyzD = ixyz[iD]
            depthB = depth[iB]
            depthC = depth[iC]
            depthD = depth[iD]
            
            # Make dictionaries to id distances
            dictB = {ixyzB[i][0][0]:depthB[i] for i in range(len(ixyzB))}
            dictC = {ixyzC[i][0][0]:depthC[i] for i in range(len(ixyzC))}
            dictD = {ixyzD[i][0][0]:depthD[i] for i in range(len(ixyzD))}
            # ID the depths
            for keyB in dictB:
                for keyC in dictC:
                    if keyB==keyC:
                        b2 = dictB[keyB]
                        c2 = dictC[keyC]
                for keyD in dictD:
                    if keyB==keyD:
                        b1 = dictB[keyB]
                        d1 = dictD[keyD]
            for keyC in dictC:
                for keyD in dictD:
                    if keyC==keyD:
                        c1 = dictC[keyC]
                        d2 = dictD[keyD]
            Bleg = b1 * (self.dx - c1) / (self.dx - b1)
            Cleg = c2 * (self.dx - d1) / (self.dx - c2)
            Dleg = d2 * (self.dx - b2) / (self.dx - d2)
            
            BigPyramid = (self.dx+Bleg) * (self.dx+Cleg) * (self.dx+Dleg) / 6
            BPyramid = Bleg * b1 * b2 / 6
            CPyramid = Cleg * c1 * c2 / 6
            DPyramid = Dleg * d1 * d2 / 6
            
            V = BigPyramid - BPyramid - CPyramid - DPyramid
            
            if V > (5/6) * self.dx**3:
                print('WARNING: Volume 4 points, 2nd case, too big')
            elif V < (1/6) * self.dx**3:
                print('WARNING: Volume 4 points, 2nd case, too small')
        return V
    
    def Points1point(self, r, point, pointneighbour):
        return self.interpolate_radius(r, point[0], pointneighbour[0])
    
    def Points2points(self, r, points, pointneighbours):  # Prism, with a right triangle base
        point1 = points[0]
        point2 = points[1]
        pointneighbour1 = [neighbour for neighbour in pointneighbours[0] if list(neighbour)!=list(point2)]
        pointneighbour2 = [neighbour for neighbour in pointneighbours[1] if list(neighbour)!=list(point1)]
        
        points = self.interpolate_radius(r, point1, pointneighbour1)
        points += self.interpolate_radius(r, point2, pointneighbour2)
        return points
    
    def Points3points(self, r, points, pointneighbours):
        # This is like the heptahedron case for the 4 point case
        # But c1 = c2 and the box width is smaller = dist(A, C)
        point1 = points[0]
        point2 = points[1]
        point3 = points[2]
        pointneighbour1 = [neighbour for neighbour in pointneighbours[0] if list(neighbour)!=list(point2) and list(neighbour)!=list(point3)]
        pointneighbour2 = [neighbour for neighbour in pointneighbours[1] if list(neighbour)!=list(point1) and list(neighbour)!=list(point3)]
        pointneighbour3 = [neighbour for neighbour in pointneighbours[2] if list(neighbour)!=list(point1) and list(neighbour)!=list(point2)]
        
        points = self.interpolate_radius(r, point1, pointneighbour1)
        points += self.interpolate_radius(r, point2, pointneighbour2)
        points += self.interpolate_radius(r, point3, pointneighbour3)
        return points
    
    def Points4points(self, r, points, pointneighbours):
        point1 = points[0]
        point2 = points[1]
        point3 = points[2]
        point4 = points[3]
        
        pointneighbour1 = [neighbour for neighbour in pointneighbours[0] if list(neighbour)!=list(point2) and list(neighbour)!=list(point3) and list(neighbour)!=list(point4)]
        pointneighbour2 = [neighbour for neighbour in pointneighbours[1] if list(neighbour)!=list(point1) and list(neighbour)!=list(point3) and list(neighbour)!=list(point4)]
        pointneighbour3 = [neighbour for neighbour in pointneighbours[2] if list(neighbour)!=list(point1) and list(neighbour)!=list(point2) and list(neighbour)!=list(point4)]
        pointneighbour4 = [neighbour for neighbour in pointneighbours[3] if list(neighbour)!=list(point1) and list(neighbour)!=list(point2) and list(neighbour)!=list(point3)]
        
        points = self.interpolate_radius(r, point1, pointneighbour1)
        points += self.interpolate_radius(r, point2, pointneighbour2)
        points += self.interpolate_radius(r, point3, pointneighbour3)
        points += self.interpolate_radius(r, point4, pointneighbour4)
        return points
    