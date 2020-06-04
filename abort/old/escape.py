

class fuel():
	def __init__(self, capsule_mass, Isp, gamma, cp, T_c, M_W, R_A=8314, g0=9.80665):
		self.R_a = R_A
		self.Isp = Isp
		self.Ve = Isp * g0
		self.gm = gamma
		self.cp = cp
		self.Tc = T_c
		self.Mw = M_W
		self.pressure_ratio()
		self.thrust(capsule_mass)
		self.mass_flow()

	def mass_flow(self):
		self.m_flow = self.F_T / self.Ve

	def pressure_ratio(self):
		# p_e / p_c
		a = (self.gm - 1) / self.gm		# (γ-1) / γ
		b = 2*self.cp					# (2γ) / (γ - 1) * Ra / Mw = 2*c_p
		self.p_e_c = (1 - (self.Ve ** 2 / (b * self.Tc))) ** (1 / a)

	def thrust(self, m, range=60, time=3):
		a = range / time ** 2
		self.F_T = m * a


hydrazine = fuel(Isp=327, gamma=1.25, cp=25.23, T_c=3250, M_W=32.0452, capsule_mass=12500)
print(hydrazine.p_e_c)
print(hydrazine.m_flow)