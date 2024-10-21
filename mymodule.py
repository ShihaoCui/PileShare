"""
This module is used to perform dispersion analysis for a cylindrical pile. 
developed by Shihao Cui, Polytechnique Montreal.
Github link: https://github.com/ShihaoCui/Blog_of_Pile.github.io
"""

@dataclass
class DPR:
    ''' 
    Object for the dispersion analysis of a cylindrical pile with a finite length L an a finite radius r

    :param fn: A mandatory (floating) vector that defines the given frequencies (Hz) which refers to the resonace frequency observed by the responses 
    :param E:    A mandatory (floating) scalar that defines the Young's modulus (Pa) of the test pile 
    :param vv:   A mandatory (floating) scalar that defines the Poisson's ratio of the test pile
    :param rho:  A mandatory (floating) scalar that defines the density (kg/m^3) for test pile
    :param r:  A mandatory (floating) scalar that defines the radius of the test pile
    :param m:  A mandatory (floating) scalar that defines the which mode would be used, where m = 0 longitudinal mode; m=1 flexural mode; m is defaluted as 0(longitudinal mode).
    :return:  Dispersion relations, i.e., the relation between the phase velocity and the given resonance frequencies
    ''' 
    
    def __init__(self, fn, E, rho, vv, r,m=0):
        self.fn = fn
        self.E = E
        self.rho = rho
        self.vv = vv
        self.r = r
        self.m = m

    def model_L(self, i): 
        ''' Define dispersion modelling and root searching algorithm for the longtudinal guided wave model for a cylinder waveguide(Pile)'''
                
        # Bessel functions at order 0, 1, and 2 of the 1st kind
        def J0(y):
            return sp.jv(0, y)
        def J1(y):
            return sp.jv(1, y)
        def J2(y):
            return sp.jv(2, y)
    
        def fun(k): 
            E = self.E
            rho  = self.rho
            vv = self.vv
            
            r = self.r
            
            cb = np.sqrt(E/rho) # bar velocity

        #     #lame coefficients
            nu = E/(2*(1+vv))
            lbd = E*vv/((1+vv)*(1-2*vv)) 

            cs = cb*np.sqrt(1/(2*(1+vv))) # shear wave speed
            cp = cb*np.sqrt((1-vv)/((1+vv)*(1-2*vv))) # longtudinal wave speed
        #     print(str(cs)+str("      ")+str(cp))

            ks = np.sqrt(np.complex(-k**2+(w**2)/(cs**2))) # the wavenumebr of the shear waves
            kp = np.sqrt(np.complex(-k**2+(w**2)/(cp**2))) # the wavenumber of the longitudinal waves

            M = np.zeros([2,2],dtype = complex)# the matrix of the disperison relation, |M|=0 is the spectrum relation
            M[0,0] = -(lbd*k**2+lbd*kp**2+2*nu*kp**2)*J0(kp*r)+2*nu*kp*J1(kp*r)/r
            M[0,1] = -2*nu*1j*k*ks*J0(ks*r)+2*nu*1j*k*J1(ks*r)/r
            M[1,0] = -2*1j*nu*k*kp*J1(kp*r)
            M[1,1] = nu*(k**2-ks**2)*J1(ks*r)
    
            return np.linalg.det(M)/1e9
    
    
        w = i*2*np.pi # w omega: angular frequency
        incre = 1000 # initialization of the phase velocity
        root = 0.00001 # the initialization of the root
        step_length = 1 # the step length of the searching algorithm
        
        for j in range(10**6):
            past = incre
            val1 = fun(w/incre)
            incre =  incre + step_length
            val2 = fun(w/incre)
            if (np.real(val1) * np.real(val2) <= 0):
                root =  optimize.ridder(fun,w/incre,w/past) # the classical Ridders' algorithm
                break 
        return (w/root)     #give one value at a frequency 
    
#     def model_F(self, i): 
#         return (w/root)     #give one value at a frequency 

    def run(self):
        ''' Define function to be used for the Parallel computing'''
        omega = np.array(self.fn) # fn is the given input frequency (Hz) vector
        root_all = [] # the result of the phase velocity by the theoretical model corresponding to the give frequency vector (fn)
        for i in range(len(fn)):
          if self.m ==0:
            root_all.append(self.model_L(omega[i]))
          elif self.m == 1:
            root_all.append(self.model_F(omega[i]))
          else:
            print("Input error; m should be either o or 1; m=0 indicating longitudinal mode,m=1 indicating flexural mode; ")
        return np.array(root_all) 