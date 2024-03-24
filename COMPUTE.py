import numpy as np
import math as math
# %% DERRIVE INFLUENCE COEFFICIENTS
# The concept of solving potential flow is rather straight forward. The linear 
# super position of source and doublet flows allows for more complex flows to
# be generated. Using the potential flow solutions for both source and doublet
# flows, the influence coefficients can be determined. 
#
# ORDER OF OPERATIONS
# [1] - Define the potential solutions for both the source and doublet elements
# [2] - Calculate for each panel the summed potential for source/doublet strength
# of 1
# [3] - The potential on a collocation by its own panel is a known constant
# [4] - Ensure that the boundaries are adhered to



# v1: [i,j] --> [ib,jb]
# v2: [i,j+1] --> [ib,jb+1]
# v3: [i+1,j+1] --> [ib+1,jb+1]
# v4: [i+1,j] --> [ib+1,jb]




# %% Initialise the parameters
m = 0.02
p = 0.4
nc = 20
ns = 10
t = 0.12
c = 1
b = 4


XB = np.zeros((2*nc-1, ns))
YB = np.zeros((2*nc-1, ns))
ZB = np.zeros((2*nc-1, ns))

XC = np.zeros((2*nc-2, ns-1))
YC = np.zeros((2*nc-2, ns-1))
ZC = np.zeros((2*nc-2, ns-1))

lx = np.zeros((2*nc-2, ns-1))
ly = np.zeros((2*nc-2, ns-1))
lz = np.zeros((2*nc-2, ns-1))

px = np.zeros((2*nc-2, ns-1))
py = np.zeros((2*nc-2, ns-1))
pz = np.zeros((2*nc-2, ns-1))

ox = np.zeros((2*nc-2, ns-1))
oy = np.zeros((2*nc-2, ns-1))
oz = np.zeros((2*nc-2, ns-1))

vx = np.ones((2*nc-2, ns-1))
vy = np.zeros((2*nc-2, ns-1))
vz = np.zeros((2*nc-2, ns-1))

ivx = np.zeros((2*nc-2, ns-1))
ivy = np.zeros((2*nc-2, ns-1))
ivz = np.zeros((2*nc-2, ns-1))

qx = np.zeros((2*nc-2, ns-1))
qy = np.zeros((2*nc-2, ns-1))
qz = np.zeros((2*nc-2, ns-1))



A = np.zeros(((2*nc-2)*(ns-1), (2*nc-2)*(ns-1)))
B = np.zeros(((2*nc-2)*(ns-1), (2*nc-2)*(ns-1)))



Sphi = np.zeros(((2*nc-2)*(ns-1), (2*nc-2)*(ns-1)))
Dphi = np.zeros(((2*nc-2)*(ns-1), (2*nc-2)*(ns-1)))





# %% Calculate the boundary points

for j in range(ns):
    x = (1-np.cos(np.linspace(0, np.pi, nc)))/2
    y = b*(1-np.cos(np.linspace(0, np.pi, ns)))/2
    
       
    # Camber definition
    Zc_1 = ((m*c) / (p**2)) * (2 * p *x[0:int(np.round(nc * p))] - x[0:int(np.round(nc * p))]**2)
    Zc_2 = ((m*c) / ((1-p)**2)) * (1-2 * p + 2 * p *x[int(np.round(nc * p)):nc] - x[int(np.round(nc * p)):nc]**2)
    camberZ = np.concatenate((Zc_1, Zc_2))
    
    # Gradient definition
    dZc_1 = ((2*m) / (p**2)) * (p - x[0:int(np.round(nc * p))])
    dZc_2 = ((2*m) / ((1-p)**2)) * (p - x[int(np.round(nc * p)):nc])
    dZc = np.concatenate((dZc_1, dZc_2))
       
    # Thickness distribution
    thicknessZ = ((c*t)/(0.2))*(0.2969*x**0.5 - 0.126*x - 0.3516*x**2 + 0.2843*x**3 - 0.1036*x**4)
    theta = np.arctan(dZc)
       
    # Evaluate the upper and lower surfaces
    Xu = x - thicknessZ*np.sin(theta)
    Xl = x + thicknessZ*np.sin(theta)
    Zu = camberZ + thicknessZ*np.cos(theta)
    Zl = camberZ - thicknessZ*np.cos(theta)
    
    # Evaluation of the boundary points
    XB[:,j] = c*np.concatenate((np.flip(Xl[1:len(Xl)]), Xu))
    ZB[:,j] = np.concatenate((np.flip(Zl[1:len(Zl)]), Zu)) 
    YB[:,j] = y[j]*np.ones(2*nc-1)


# %% Calculate the collocation points and panel areas
for j in range(ns-1):    
    for i in range(2*nc-2):

        
# v1: [i,j] --> [ib,jb]
# v2: [i,j+1] --> [ib,jb+1]
# v3: [i+1,j+1] --> [ib+1,jb+1]
# v4: [i+1,j] --> [ib+1,jb]


        # Calculate the collocation points
        XC[i,j] = (XB[i,j] + XB[i+1,j] + XB[i,j+1] + XB[i+1,j+1])/4
        YC[i,j] = (YB[i,j] + YB[i+1,j] + YB[i,j+1] + YB[i+1,j+1])/4
        ZC[i,j] = (ZB[i,j] + ZB[i+1,j] + ZB[i,j+1] + ZB[i+1,j+1])/4
        
        # Calculate the x-wise normal vector
        lx[i,j] = (XB[i,j] + XB[i,j+1] - XB[i+1,j+1] - XB[i+1,j])/2
        ly[i,j] = (YB[i,j] + YB[i,j+1] - YB[i+1,j+1] - YB[i+1,j])/2
        lz[i,j] = (ZB[i,j] + ZB[i,j+1] - ZB[i+1,j+1] - ZB[i+1,j])/2
        
        # Calculate the y-wise normal vector
        px[i,j] = (XB[i,j+1] + XB[i+1,j+1] - XB[i+1,j] - XB[i,j])/2
        py[i,j] = (YB[i,j+1] + YB[i+1,j+1] - YB[i+1,j] - YB[i,j])/2
        pz[i,j] = (ZB[i,j+1] + ZB[i+1,j+1] - ZB[i+1,j] - ZB[i,j])/2
        
        # Calculate the z-wise normal vector
        d1 = np.array((XB[i+1,j+1] - XB[i,j], YB[i+1,j+1] - YB[i,j], ZB[i+1,j+1] - ZB[i,j]))
        d2 = np.array((XB[i,j+1] - XB[i+1,j], YB[i,j+1] - YB[i+1,j], ZB[i,j+1] - ZB[i+1,j]))
        d3 = np.cross(d1,d2)
        
        
        ox[i,j] = d3[0]
        oy[i,j] = d3[1]
        oz[i,j] = d3[2]

# %% Calculate the total number of panels

numPans = len(XC[0,:])*len(XC[:,0])

# %% ELEMENTS FOR PROCESSING
#
#

def ELEMENTS(XB, YB, ZB, XC, YC, ZC, ib, jb):

    d12 = np.sqrt((XB[ib,jb+1] - XB[ib,jb])**2 + (YB[ib,jb+1] - YB[ib,jb])**2) 
    d23 = np.sqrt((XB[ib+1,jb+1] - XB[ib,jb+1])**2 + (YB[ib+1,jb+1] - YB[ib,jb+1])**2) 
    d34 = np.sqrt((XB[ib+1,jb] - XB[ib+1,jb+1])**2 + (YB[ib+1,jb] - YB[ib+1,jb+1])**2)
    d41 = np.sqrt((XB[ib,jb] - XB[ib+1,jb])**2 + (YB[ib,jb] - YB[ib+1,jb])**2) 
    
    if (XB[ib,jb+1] - XB[ib,jb]) <= 0:
        m12 = 0
    else:
        m12 = (YB[ib,jb+1] - YB[ib,jb])/(XB[ib,jb+1] - XB[ib,jb])

    m23 = (YB[ib+1,jb+1] - YB[ib,jb+1])/(XB[ib+1,jb+1] - XB[ib,jb+1])

    if (XB[ib+1,jb] - XB[ib+1,jb+1]) <= 0:
        m34 = 0
    else:
        m34 = (YB[ib+1,jb] - YB[ib+1,jb+1])/(XB[ib+1,jb] - XB[ib+1,jb+1])

    m41 = (YB[ib,jb] - YB[ib+1,jb])/(XB[ib,jb] - XB[ib+1,jb])
    
    r1 = np.sqrt((XC[ic,jc] - XB[ib,jb])**2 + (YC[ic,jc] - YB[ib,jb])**2 + ZC[ic,jc]**2)
    r2 = np.sqrt((XC[ic,jc] - XB[ib,jb+1])**2 + (YC[ic,jc] - YB[ib,jb+1])**2 + ZC[ic,jc]**2)
    r3 = np.sqrt((XC[ic,jc] - XB[ib+1,jb+1])**2 + (YC[ic,jc] - YB[ib+1,jb+1])**2 + ZC[ic,jc]**2)
    r4 = np.sqrt((XC[ic,jc] - XB[ib+1,jb])**2 + (YC[ic,jc] - YB[ib+1,jb])**2 + ZC[ic,jc]**2)
    
    e1 = (XC[ic,jc] - XB[ib,jb])**2 + ZC[ic,jc]**2
    e2 = (XC[ic,jc] - XB[ib,jb+1])**2 + ZC[ic,jc]**2
    e3 = (XC[ic,jc] - XB[ib+1,jb+1])**2 + ZC[ic,jc]**2
    e4 = (XC[ic,jc] - XB[ib+1,jb])**2 + ZC[ic,jc]**2
    
    h1 = (XC[ic,jc] - XB[ib,jb])*(YC[ic,jc] - YB[ib,jb])  
    h2 = (XC[ic,jc] - XB[ib,jb+1])*(YC[ic,jc] - YB[ib,jb+1]) 
    h3 = (XC[ic,jc] - XB[ib+1,jb+1])*(YC[ic,jc] - YB[ib+1,jb+1]) 
    h4 = (XC[ic,jc] - XB[ib+1,jb])*(YC[ic,jc] - YB[ib+1,jb]) 
    
    return d12, d23, d34, d41, m12, m23, m34, m41, r1, r2, r3, r4, e1, e2, e3, e4, h1, h2, h3, h4



# %% QUADRILATERAL SOURCE POTENTIAL
# Calculates: Influence of the source potential from the boundary points (panels) on the collocation points (measurements)
# Uses: Hess and Smith derrivation for the planar element
# Produces: A matrix influence coefficients

for ic in range(len(XC[:,0])):
    for jc in range(len(XC[0,:])):
        for ib in range(len(XC[:,0])):
            for jb in range(len(XC[0,:])):
            
                d12, d23, d34, d41, m12, m23, m34, m41, r1, r2, r3, r4, e1, e2, e3, e4, h1, h2, h3, h4 = ELEMENTS(XB, YB, ZB, XC, YC, ZC, ib, jb)

                Sterm1 = ((XC[ic,jc]-XB[ib,jb])*(YB[ib,jb+1]-YB[ib,jb]) - (YC[ic,jc]-YB[ib,jb])*(XB[ib,jb+1]-XB[ib,jb]))*(np.log((r1+r2+d12)/(r1+r2-d12)))
                Sterm2 = ((XC[ic,jc]-XB[ib,jb+1])*(YB[ib+1,jb+1]-YB[ib,jb+1]) - (YC[ic,jc]-YB[ib,jb+1])*(XB[ib+1,jb+1]-XB[ib,jb+1]))*(np.log((r2+r3+d23)/(r2+r3-d23)))
                Sterm3 = ((XC[ic,jc]-XB[ib+1,jb+1])*(YB[ib+1,jb]-YB[ib+1,jb+1]) - (YC[ic,jc]-YB[ib+1,jb+1])*(XB[ib+1,jb]-XB[ib+1,jb+1]))*(np.log((r3+r4+d34)/(r3+r4-d34)))
                Sterm4 = ((XC[ic,jc]-XB[ib+1,jb])*(YB[ib,jb]-YB[ib+1,jb]) - (YC[ic,jc]-YB[ib+1,jb])*(XB[ib,jb]-XB[ib+1,jb]))*(np.log((r4+r1+d41)/(r4+r1-d41)))
                
                Sterm5 = math.atan2((m12*e1 - h1),(ZC[ic,jc]*r1)) - math.atan2((m12*e2 - h2),(ZC[ic,jc]*r2))
                Sterm6 = math.atan2((m23*e2 - h2),(ZC[ic,jc]*r2)) - math.atan2((m23*e3 - h2),(ZC[ic,jc]*r3))
                Sterm7 = math.atan2((m34*e3 - h3),(ZC[ic,jc]*r3)) - math.atan2((m34*e4 - h2),(ZC[ic,jc]*r4))
                Sterm8 = math.atan2((m41*e4 - h4),(ZC[ic,jc]*r4)) - math.atan2((m41*e1 - h2),(ZC[ic,jc]*r1))
                
                # # Account for influence of own measured panel
                # if jc*len(XC[:,0]) + ic == jb*len(XC[:,0]) + ib:
                #     if oz[ic,jc] > 0:
                #         Sphi[jb*len(XC[:,0]) + ib, jc*len(XC[:,0]) + ic] = 0.5
                #     else:
                #         Sphi[jb*len(XC[:,0]) + ib, jc*len(XC[:,0]) + ic] = -0.5
                # else:
                    
                Sphi[jb*len(XC[:,0]) + ib, jc*len(XC[:,0]) + ic] = (-1)/(4*np.pi) *((Sterm1 + Sterm2 + Sterm3 + Sterm4) - abs(ZC[ic,jc])*(Sterm5 + Sterm6 + Sterm7 + Sterm8))
        
# %% QUADRILATERAL DOUBLET POTENTIAL
# Calculates: Influence of the doublet potential from the boundary points (panels) on the collocation points (measurements)
# Uses: Hess and Smith derrivation for the planar element
# Produces: A matrix influence coefficients 
             
for ic in range(len(XC[:,0])):
    for jc in range(len(XC[0,:])):
        for ib in range(len(XC[:,0])):
            for jb in range(len(XC[0,:])):
                    
                d12, d23, d34, d41, m12, m23, m34, m41, r1, r2, r3, r4, e1, e2, e3, e4, h1, h2, h3, h4 = ELEMENTS(XB, YB, ZB, XC, YC, ZC, ib, jb)

                Dterm1 = math.atan2((m12*e1 - h1),(ZC[ic,jc]*r1)) - math.atan2((m12*e2 - h2),(ZC[ic,jc]*r2))
                Dterm2 = math.atan2((m23*e2 - h2),(ZC[ic,jc]*r2)) - math.atan2((m23*e3 - h2),(ZC[ic,jc]*r3))
                Dterm3 = math.atan2((m34*e3 - h3),(ZC[ic,jc]*r3)) - math.atan2((m34*e4 - h2),(ZC[ic,jc]*r4))
                Dterm4 = math.atan2((m41*e4 - h4),(ZC[ic,jc]*r4)) - math.atan2((m41*e1 - h2),(ZC[ic,jc]*r1))
                
                
                Dphi[jb*len(XC[:,0]) + ib, jc*len(XC[:,0]) + ic] = (1)/(4*np.pi) *(Dterm1 + Dterm2 + Dterm3 + Dterm4)
                
                # Apply the 3D Kutta condition for the fictious wake panel
                if ic == 0:
                    DphiW1 = Dphi[jb*len(XC[:,0]) + ib, jc*len(XC[:,0]) + ic]
                elif ic == len(XC[:,0])-1:
                    DphiW2 = Dphi[jb*len(XC[:,0]) + ib, jc*len(XC[:,0]) + ic]
                    DphiW = DphiW2 - DphiW1
                    
                    # Apply to only the trailing edge panels
                    Dphi[jb*len(XC[:,0]) + ib, jc*len(XC[:,0]) + ic] = Dphi[jb*len(XC[:,0]) + ib, jc*len(XC[:,0]) + ic] + DphiW 
                    Dphi[jb*len(XC[:,0]) + ib, jc*len(XC[:,0])] = Dphi[jb*len(XC[:,0]) + ib, jc*len(XC[:,0])] - DphiW 
                
                
# %% SOLVE LINEAR SET OF EQUATIONS
#
#


# Dirichlet Boundary Condition for the Source Strengths
sigma = np.dot(np.array((np.ravel(ox), np.ravel(oy), np.ravel(oz))).T, 
               np.array((np.ravel(vx), np.ravel(vy), np.ravel(vz))))[:,0]














# Populate the A matrix
for i in range(numPans):                                                       
    for j in range(numPans):            
        if (i == j):                                                           
            A[i,j] = 0.5   
        else:                                                                 
            A[i,j] = Dphi[i,j]    

B = np.dot(sigma, Sphi)


mu = np.linalg.solve(A,B)



#mu = np.reshape(mu, (2*nc-2, ns-1))
# #sigma = np.reshape(resArrSig, (2*nc-2, ns-1))

mu = np.rot90(np.reshape(mu, (ns-1,2*nc-2)))
sigma = np.reshape(sigma, (2*nc-2, ns-1))


# Calculate the singularity flows
for jb in range(len(XC[0,:])-1):
    for ib in range(len(XC[:,0])-1):
        qx[ib,jb] = - ((mu[ib+1,jb] - mu[ib,jb])/(lx[ib+1,jb] - lx[ib,jb]))
        qy[ib,jb] = - ((mu[ib,jb+1] - mu[ib,jb])/(py[ib,jb+1] - py[ib,jb]))
        
qz = -sigma

# calculate the net induced velocities
ivx = vx + qx
ivy = vy + qy
ivz = vz + qz














# # Populate the B matrix
# for i in range(numPans):
    



#  # Populate the B matrix
#  for i in range(numPans):
#      B[i, k] = -Vinf*2*np.pi*np.cos(beta[k,i])
 
 
#  # Include the Kutta Condition on the bottom of the A matrix
#  for i in range(foils):
#      for j in range(numPans):
#          A[numPans+i, j, k] = J[startInd[i],j, k] + J[endInd[i]-foils,j, k]
#      A[numPans+i,numPans+i, k] = -sum(L[startInd[i],:, k] + L[endInd[i]-foils,:, k]) + 2*np.pi
     
#  # Include the Kutta Condition at the final point of the B matrix
#  for i in range(foils):
#      B[numPans+i, k] = -Vinf*2*np.pi*(np.sin(beta[k, startInd[i]]) + np.sin(beta[k, endInd[i]-foils]))





# resArrSig = np.ravel(ox) * np.ravel(vx) + np.ravel(oy) * np.ravel(vy) + np.ravel(oz) * np.ravel(vz)
# B = -np.dot(Sphi,resArrSig)

# # %% BUILD THE LHS
# #
# #

# A = Dphi

# # %% SOLVE
# #
# #

# resArrMu = np.linalg.solve(A, B)

# # Reshape into useful form
# mu = np.reshape(resArrMu, (2*nc-2, ns-1))
# #sigma = np.reshape(resArrSig, (2*nc-2, ns-1))

# #mu = np.rot90(np.reshape(resArrMu, (ns-1,2*nc-2)))
# sigma = np.reshape(resArrSig, (2*nc-2, ns-1))

# # %% INDUCED TOTAL VELOCITY
# #
# #


# # Calculate the singularity flows
# for jb in range(len(XC[0,:])-1):
#     for ib in range(len(XC[:,0])-1):
#         qx[ib,jb] = - ((mu[ib+1,jb] - mu[ib,jb])/(lx[ib+1,jb] - lx[ib,jb]))
#         qy[ib,jb] = - ((mu[ib,jb+1] - mu[ib,jb])/(py[ib,jb+1] - py[ib,jb]))
        
# qz = -sigma

# # calculate the net induced velocities
# ivx = vx + qx
# ivy = vy + qy
# ivz = vz + qz

# cp = 1 - (sum(np.array((ivx, ivy, ivz)))**2)/(sum(np.array((vx, vy, vz)))**2)






























# %% INDUCED DOUBLET VELOCITIES
#
#






# # Calculate doublet induced x velocity: ivx

# for jb in range(len(XC[0,:])):
#     for ib in range(len(XC[:,0])):
    
#         d12, d23, d34, d41, m12, m23, m34, m41, r1, r2, r3, r4, e1, e2, e3, e4, h1, h2, h3, h4 = ELEMENTS(XB, YB, ZB, XC, YC, ZC, ib, jb)


#         Sterm1 = ZC[ib,jb]*((YB[ib,jb] - YB[ib,jb+1])*(r1+r2))
#         Sterm2 = r1*r2*(r1*r2 - ((XC[ib,jb] - XB[ib,jb])*(XC[ib,jb] - XB[ib,jb+1]) + (YC[ib,jb] - YB[ib,jb])*(YC[ib,jb] - YB[ib,jb+1]) + (ZC[ib,jb]**2)))
        
#         Sterm3 = ZC[ib,jb]*((YB[ib,jb+1] - YB[ib+1,jb+1])*(r2+r3))
#         Sterm4 = r2*r3*(r2*r3 - ((XC[ib,jb] - XB[ib,jb+1])*(XC[ib,jb] - XB[ib+1,jb+1]) + (YC[ib,jb] - YB[ib,jb+1])*(YC[ib,jb] - YB[ib+1,jb+1]) + (ZC[ib,jb]**2)))
        
#         Sterm5 = ZC[ib,jb]*((YB[ib+1,jb+1] - YB[ib+1,jb])*(r3+r4))
#         Sterm6 = r3*r4*(r3*r4 - ((XC[ib,jb] - XB[ib+1,jb+1])*(XC[ib,jb] - XB[ib+1,jb]) + (YC[ib,jb] - YB[ib+1,jb+1])*(YC[ib,jb] - YB[ib+1,jb]) + (ZC[ib,jb]**2)))
        
#         Sterm7 = ZC[ib,jb]*((YB[ib+1,jb] - YB[ib,jb])*(r4+r1))
#         Sterm8 = r4*r1*(r4*r1 - ((XC[ib,jb] - XB[ib+1,jb])*(XC[ib,jb] - XB[ib,jb]) + (YC[ib,jb] - YB[ib+1,jb])*(YC[ib,jb] - YB[ib,jb]) + (ZC[ib,jb]**2)))
        
#         ivx[ib,jb] = (mu[ib,jb]/(4*np.pi))*(Sterm1/Sterm2 + Sterm3/Sterm4 + Sterm5/Sterm6 + Sterm7/Sterm8)

# # Calculate doublet induced y velocity: ivy
# for jc in range(len(XC[0,:])):
#     for ic in range(len(XC[:,0])):
#         for jb in range(len(XC[0,:])):
#             for ib in range(len(XC[:,0])):
            
#                 d12, d23, d34, d41, m12, m23, m34, m41, r1, r2, r3, r4, e1, e2, e3, e4, h1, h2, h3, h4 = ELEMENTS(XB, YB, ZB, XC, YC, ZC, ib, jb)


#                 Sterm1 = ZC[ic,jc]*((YB[ib,jb+1] - YB[ib,jb])*(r1+r2))
#                 Sterm2 = r1*r2*(r1*r2 - ((XC[ic,jc] - XB[ib,jb])*(XC[ib,jb] - XB[ib,jb+1]) + (YC[ic,jc] - YB[ib,jb])*(YC[ib,jb] - YB[ib,jb+1]) + (ZC[ic,jc]**2)))
                
#                 Sterm3 = ZC[ic,jc]*((YB[ib+1,jb+1] - YB[ib,jb+1])*(r2+r3))
#                 Sterm4 = r2*r3*(r2*r3 - ((XC[ic,jc] - XB[ib,jb+1])*(XC[ic,jc] - XB[ib+1,jb+1]) + (YC[ic,jc] - YB[ib,jb+1])*(YC[ic,jc] - YB[ib+1,jb+1]) + (ZC[ic,jc]**2)))
                
#                 Sterm5 = ZC[ic,jc]*((YB[ib+1,jb] - YB[ib+1,jb+1])*(r3+r4))
#                 Sterm6 = r3*r4*(r3*r4 - ((XC[ic,jc] - XB[ib+1,jb+1])*(XC[ic,jc] - XB[ib+1,jb]) + (YC[ic,jc] - YB[ib+1,jb+1])*(YC[ic,jc] - YB[ib+1,jb]) + (ZC[ic,jc]**2)))
                
#                 Sterm7 = ZC[ic,jc]*((YB[ib,jb] - YB[ib+1,jb])*(r4+r1))
#                 Sterm8 = r4*r1*(r4*r1 - ((XC[ic,jc] - XB[ib+1,jb])*(XC[ic,jc] - XB[ib,jb]) + (YC[ic,jc] - YB[ib+1,jb])*(YC[ic,jc] - YB[ib,jb]) + (ZC[ic,jc]**2)))
                
#         ivy[ic,jc] = (mu[ic,jc]/(4*np.pi))*(Sterm1/Sterm2 + Sterm3/Sterm4 + Sterm5/Sterm6 + Sterm7/Sterm8)              

# # Calculate doublet induced z velocity: ivz
# for jc in range(len(XC[0,:])):
#     for ic in range(len(XC[:,0])):
#         for jb in range(len(XC[0,:])):
#             for ib in range(len(XC[:,0])):
            
#                 d12, d23, d34, d41, m12, m23, m34, m41, r1, r2, r3, r4, e1, e2, e3, e4, h1, h2, h3, h4 = ELEMENTS(XB, YB, ZB, XC, YC, ZC, ib, jb)


#                 Sterm1 = ((XC[ic,jc] - XB[ib,jb+1])*(YC[ic,jc] - YB[ib,jb]) - (XC[ic,jc] - XB[ib,jb])*(YC[ic,jc] - YB[ib,jb+1]))*(r1+r2)
#                 Sterm2 = r1*r2*(r1*r2 - ((XC[ic,jc] - XB[ib,jb])*(XC[ib,jb] - XB[ib,jb+1]) + (YC[ic,jc] - YB[ib,jb])*(YC[ib,jb] - YB[ib,jb+1])))
               
#                 Sterm3 = ((XC[ic,jc] - XB[ib+1,jb+1])*(YC[ic,jc] - YB[ib,jb+1]) - (XC[ic,jc] - XB[ib,jb+1])*(YC[ic,jc] - YB[ib+1,jb+1]))*(r2+r3)
#                 Sterm4 = r2*r3*(r2*r3 - ((XC[ic,jc] - XB[ib,jb+1])*(XC[ic,jc] - XB[ib+1,jb+1]) + (YC[ic,jc] - YB[ib,jb+1])*(YC[ic,jc] - YB[ib+1,jb+1])))
            
#                 Sterm5 = ((XC[ic,jc] - XB[ib+1,jb])*(YC[ic,jc] - YB[ib+1,jb+1]) - (XC[ic,jc] - XB[ib+1,jb+1])*(YC[ic,jc] - YB[ib+1,jb]))*(r3+r4)
#                 Sterm6 = r3*r4*(r3*r4 - ((XC[ic,jc] - XB[ib+1,jb+1])*(XC[ic,jc] - XB[ib+1,jb]) + (YC[ic,jc] - YB[ib+1,jb+1])*(YC[ic,jc] - YB[ib+1,jb])))
   
#                 Sterm7 = ((XC[ic,jc] - XB[ib,jb])*(YC[ic,jc] - YB[ib+1,jb]) - (XC[ic,jc] - XB[ib+1,jb])*(YC[ic,jc] - YB[ib,jb]))*(r4+r1)
#                 Sterm8 = r4*r1*(r4*r1 - ((XC[ic,jc] - XB[ib+1,jb])*(XC[ic,jc] - XB[ib,jb]) + (YC[ic,jc] - YB[ib+1,jb])*(YC[ic,jc] - YB[ib,jb])))
                
#         ivz[ic,jc] = (mu[ic,jc]/(4*np.pi))*(Sterm1/Sterm2 + Sterm3/Sterm4 + Sterm5/Sterm6 + Sterm7/Sterm8)                          




























































































        
        
        
        # Sterm5 = np.arctan((m12*e1 - h1)/(ZC[ib,jb]*r1)) - np.arctan((m12*e2 - h2)/(ZC[ib,jb]*r2))
        # Sterm6 = np.arctan((m23*e2 - h2)/(ZC[ib,jb]*r2)) - np.arctan((m23*e3 - h3)/(ZC[ib,jb]*r3))
        # Sterm7 = np.arctan((m34*e3 - h3)/(ZC[ib,jb]*r3)) - np.arctan((m34*e4 - h4)/(ZC[ib,jb]*r4))  
        # Sterm8 = np.arctan((m41*e4 - h4)/(ZC[ib,jb]*r4)) - np.arctan((m41*e1 - h1)/(ZC[ib,jb]*r1))

        # Dterm1 = np.arctan((m12*e1 - h1)/(ZC[ib,jb]*r1)) - np.arctan((m12*e2 - h2)/(ZC[ib,jb]*r2))
        # Dterm2 = np.arctan((m23*e2 - h2)/(ZC[ib,jb]*r2)) - np.arctan((m23*e3 - h3)/(ZC[ib,jb]*r3))
        # Dterm3 = np.arctan((m34*e3 - h3)/(ZC[ib,jb]*r3)) - np.arctan((m34*e4 - h4)/(ZC[ib,jb]*r4))  
        # Dterm4 = np.arctan((m41*e4 - h4)/(ZC[ib,jb]*r4)) - np.arctan((m41*e1 - h1)/(ZC[ib,jb]*r1))
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        