class Atom(object):
    def __init__(self,atomic_number,molecular_mass):
        self.atomic_number = atomic_number
        self.molecular_mass = molecular_mass
        self.ground_state_energy_level = 1
    def energize(self,energy):
        print "Current energy state is %i" %self.ground_state_energy_level
        self.ground_state_energy_level = self.ground_state_energy_level + energy
        print "Successfully energized the energy level to %s" %self.ground_state_energy_level
        

Carbon=Atom(12,"C") 
dir(Carbon)
print Carbon.atomic_number      
print Carbon.molecular_mass
print Carbon.ground_state_energy_level
Carbon.energize(2)
print Carbon.ground_state_energy_level

