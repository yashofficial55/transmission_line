import matplotlib.pyplot as plt
import numpy as np
plt.style.use("ggplot")

#-------------------------------------------------------------------------------
# Auxiliary functions

# Complex power
def complex_power(p, q):
    return np.sqrt(p**2 + q**2)
    
# Current
def current(s, V):
    return s/(np.sqrt(3)*V)

# Voltage
def voltage(s, I):
    return s/(np.sqrt(3)*I)

# Power (3 phase)
def power_3f(a, I):
    return 3 * a * I**2

# Vectorized version of round
round_ = np.vectorize(round)

#-------------------------------------------------------------------------------
# 3 phase electric load

class electric_load(object):

    def __init__(self, pf, Vn, S = None, P = None):
        # Electric load
        #
        # Power factor pf
        # Nominal voltage Vn [V]
        # Apparent (or complex) power [VA]
        # Power [W]
        #
        if S is None and P is None:
            raise ValueError("Either S or P should be assigned a value...")
        elif S is None:
            self.S = P/pf
            self.P = P
        else:
            self.S = S
            self.P = S*pf

        self.pf = pf
        self.Q = self.P * np.tan(np.arccos(pf))
        self.Vn = Vn
    
    def current_draw(self):
        # Get the current draw of the electric load
        # Output: I
        #
        return self.S/(np.sqrt(3)*self.Vn)
    
    def get_electric_load(self):
        # Return the load characteristics
        #        
        # Output: V,I,p,q
        #
        p_ = self.P
        q_ = self.Q        
        V = self.Vn
        I = self.current_draw()
        return V, I, p_, q_
    
    def add_load_in_parallel(self, p, q, V):
        # Add a load in parallel in a given line.
        # Get V,I,p,q at the top of the transmission line
        #
        # Input: p,q,V at the bottom of the line before the load
        # Output: V,I,p,q at the top of the line
        #        
        p_ = p + self.P
        q_ = q + self.Q
        s_ = complex_power(p_, q_)
        I = current(s_, V)
        return V, I, p_, q_

#-------------------------------------------------------------------------------
# Transformer class

class transformer(object):
    
    def __init__(self, V1n, V2n, Sn, vcc_pc):
        # Transfomer
        #
        # V1n = nominal voltage at primary [V]
        # V2n = nominal voltage at secondary [V]
        # Sn = nominal apparent power [VA]
        # vcc_pc = short circuit voltage [%]
        #
        self.V1n = V1n
        self.V2n = V2n
        self.Sn = Sn
        self.k = V1n/V2n
        self.vcc_pc = vcc_pc
        self.Rcc = 0
        self.Xcc = (vcc_pc * V2n**2)/(100 * Sn)
        
    def primary_values(self, p, q, I2):
        # Get the values at the primary side of the transfomer
        #
        # Input: p,q,I at the secondary side of the transfomer
        # Output: V,I,p,q at the primary side of the transfomer
        #        
        p_ = p + power_3f(self.Rcc, I2)
        q_ = q + power_3f(self.Xcc, I2)
        s_ = complex_power(p_, q_)
        V2 = voltage(s_, I2)
        V1 = self.k * V2
        I1 = I2 / self.k
        return V1, I1, p_, q_

#-------------------------------------------------------------------------------
# Transmission line class
class transmission_line(object):

    def __init__(self, R, X, L):
        # Transmission line
        # 
        # R = resistance [Ohm/km]
        # X = reactance [Ohm/km]
        # L = length [km]
        #
        self.L = L
        self.R = R * L
        self.X = X * L
        self.Z = np.sqrt(R**2 + X**2)
        #self.Z_cmp = R + 1j * X

    
    def same_current(self, p, q, I):
        # Get V,I,p,q at the top of the transmission line
        #
        # Input: p,q,I at the bottom of the line
        # Output: V,I,p,q at the top of the line
        #
        p_ = p + power_3f(self.R, I)
        q_ = q + power_3f(self.X, I)
        s_ = complex_power(p_, q_)
        V = voltage(s_, I)
        return V, I, p_, q_

#-------------------------------------------------------------------------------
# Example of application

# Loads
c1 = electric_load(pf=0.90, Vn=15000, P=1000000)
c2 = electric_load(pf=0.90, Vn=15000, S=1000000)
c3 = electric_load(pf=0.95, Vn=15000, S=2000000)

# Transmission lines
tl_2 = transmission_line(R=0.1, X=0.1, L=10)
tl_1 = transmission_line(R=0.1, X=0.1, L=5)

# Transformers
T_1 = transformer(V1n=132000, V2n=15000, Sn=25000000, vcc_pc=12.5)

# Solve transmission line
Va, Ia, Pa, Qa = c1.get_electric_load()
Vb, Ib, Pb, Qb = tl_2.same_current(Pa, Qa, Ia)
Vb, Ib, Pb, Qb = c2.add_load_in_parallel(Pb, Qb, Vb)
Vb, Ib, Pb, Qb = c3.add_load_in_parallel(Pb, Qb, Vb)
Vc, Ic, Pc, Qc = tl_1.same_current(Pb, Qb, Ib)
Vd, Id, Pd, Qd = T_1.primary_values(Pc, Qc, Ic)

# Get ready for printing out results
sections = ["a", "b", "c", "d"]
V = round_([Va, Vb, Vc, Vd], 2)
I = round_([Ia, Ib, Ic, Id], 2)

dV_ = ["Vdc", "Vcb", "Vba"]
dV = [Vd/132000 - Vc/15000, Vc/15000 - Vb/15000, Vb/15000 - Va/15000]
dV = [i * 100 for i in dV]

# Print results
print("Section", "\tVoltage [kV]", "\tCurrent [A]")
for i,k in enumerate(sections):
    print(str(k), "\t\t" + str(round(V[i]/1000, 2)), "\t\t" + str(I[i]))

print("\nVoltage drop:")
for i,k in enumerate(dV_):
    print(str(k), "\t" + str(round(dV[i], 2)) + "%")

#-------------------------------------------------------------------------------
# Plot
plt.figure()
plt.subplot(3,1,1)
plt.plot(I)
plt.title("Current [A]")
plt.xticks(range(len(sections)), sections, size='small')
plt.subplot(3,1,2)
plt.plot(V/1000)
plt.xticks(range(len(sections)), sections, size='small')
plt.title("Voltage [kV]")
plt.subplot(3,1,3)
plt.plot(dV)
plt.xticks(range(len(dV_)), dV_, size='small')
plt.title("Voltage drop %")
plt.show()
