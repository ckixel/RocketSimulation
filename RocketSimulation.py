
from math import pi, cos, sin, atan2
from __future__ import annotations
import sys
import matplotlib.pyplot as plt
import numpy as np

xpoints = np.array([])
ypoints_altitude = np.array([])
ypoints_speedX = np.array([])
ypoints_speedY = np.array([])
ypoints_localSpeedX = np.array([])
ypoints_localSpeedY = np.array([])
ypoints_fuel = np.array([])
xpoints_position = np.array([])
ypoints_position = np.array([])
zpoints_position = np.array([])

airDensitySpreadsheet = [
[0, 1.2250],
[500, 1.1673],
[1000, 1.1117],
[1500, 1.0581],
[2000, 1.0065],
[2500, 0.9569],
[3000, 0.9093],
[4000, 0.8194],
[5000, 0.7365],
[6000, 0.6601],
[7000, 0.59],
[8000, 0.5258],
[9000, 0.4671],
[10_000, 0.4135],
[11_000, 0.3648],
[12_000, 0.3119],
[14_000, 0.2279],
[16_000, 0.1665],
[18_000, 0.1216],
[20_000, 0.0889],
[24_000, 0.0469],
[28_000, 0.0251],
[32_000, 0.0136],
[36_000, 7.26 * 10 ** (-3)],
[40_000, 4.00 * 10 ** (-3)],
[50_000, 1.03 * 10 ** (-3)],
[60_000, 3.00 * 10 ** (-4)],
[80_000, 1.85 * 10 ** (-5)],
[100_000, 5.55 * 10 ** (-7)],
[150_000, 2.00 * 10 ** (-9)],
[200_000, 2.52 * 10 ** (-10)],
[300_000, 1.92 * 10 ** (-11)],
[500_000, 5.21 * 10 ** (-13)],
[700_000, 3.07 * 10 ** (-14)],
[1_000_000, 3.56 * 10 ** (-15)]
]

def getAirDensity(altitude):
	if altitude <= 0: return airDensitySpreadsheet[0][1]
	if altitude >= airDensitySpreadsheet[-1][0]: return airDensitySpreadsheet[-1][1]
	for i, el in enumerate(airDensitySpreadsheet):
		if altitude == el[0]: return el[1]
		if altitude < el[0]:
			p = airDensitySpreadsheet[i-1]
			k = (altitude - p[0]) / (el[0] - p[0])
			return p[1] * (1 - k) + el[1] * k

angleSpreadsheet = [
[5, 0],
[20, 10],
[60, 15],
[150, 20],
[160, 5],
[300, 10],
[400, 25],
[540, 40],
[550, 55],
[800, 65],
[999, 70],
[1000, 90]
]

def getAngle(time):
	if time <= 0: return angleSpreadsheet[0][1]
	if time >= angleSpreadsheet[-1][0]: return angleSpreadsheet[-1][1]
	for i, el in enumerate(angleSpreadsheet):
		if time == el[0]: return el[1]
		if time < el[0]:
			p = angleSpreadsheet[i-1]
			k = (time - p[0]) / (el[0] - p[0])
			return p[1] * (1 - k) + el[1] * k

thrustSpreadsheet = [
[0, 1],
[600, 1],
[601, 0.2],
[2000, 0.2],
[2001, 0]
]

def getThrust(time):
	if time <= 0: return thrustSpreadsheet[0][1]
	if time >= thrustSpreadsheet[-1][0]: return thrustSpreadsheet[-1][1]
	for i, el in enumerate(thrustSpreadsheet):
		if time == el[0]: return el[1]
		if time < el[0]:
			p = thrustSpreadsheet[i-1]
			k = (time - p[0]) / (el[0] - p[0])
			return p[1] * (1 - k) + el[1] * k

G = 6.6743 * 10 ** (-11)

class Vec2:
	def __init__(self, x: float, y: float):
		self.x = x
		self.y = y
	
	def __str__(self):
		return f'x: {self.x}, y: {self.y}'

	def __add__(self, other):
		x = other.x if type(other) == Vec2 else other
		y = other.y if type(other) == Vec2 else other
		return Vec2(self.x + x, self.y + y)

	def __iadd__(self, other):
		x = other.x if type(other) == Vec2 else other
		y = other.y if type(other) == Vec2 else other
		self.x += x
		self.y += y
		return self

	def __sub__(self, other):
		x = other.x if type(other) == Vec2 else other
		y = other.y if type(other) == Vec2 else other
		return Vec2(self.x - x, self.y - y)

	def __isub__(self, other):
		x = other.x if type(other) == Vec2 else other
		y = other.y if type(other) == Vec2 else other
		self.x -= x
		self.y -= y
		return self

	def __mul__(self, other):
		x = other.x if type(other) == Vec2 else other
		y = other.y if type(other) == Vec2 else other
		return Vec2(self.x * x, self.y * y)

	def __imul__(self, other):
		x = other.x if type(other) == Vec2 else other
		y = other.y if type(other) == Vec2 else other
		self.x *= x
		self.y *= y
		return self

	def __truediv__(self, other):
		x = other.x if type(other) == Vec2 else other
		y = other.y if type(other) == Vec2 else other
		if x == 0 or y == 0:
			return Vec2.zero()
		return Vec2(self.x / x, self.y / y)

	def __itruediv__(self, other):
		x = other.x if type(other) == Vec2 else other
		y = other.y if type(other) == Vec2 else other
		if x == 0 or y == 0:
			self = Vec2.zero()
			return
		self.x /= x
		self.y /= y
		return self

	def __pow__(self, other):
		x = other.x if type(other) == Vec2 else other
		y = other.y if type(other) == Vec2 else other
		return Vec2(self.x ** x, self.y ** y)

	def __ipow__(self, other):
		x = other.x if type(other) == Vec2 else other
		y = other.y if type(other) == Vec2 else other
		self.x **= x
		self.y **= y
		return self
	
	def length(self) -> float:
		return (self.x ** 2 + self.y ** 2) ** .5
	
	def normalized(self) -> Vec2:
		d = self.length()
		if d == 0: return Vec2.zero()
		return Vec2(self.x / d, self.y / d)

	def negative(self) -> Vec2:
		return Vec2(-self.x, -self.y)
	
	def rotated(self, angle: float) -> Vec2:
		return Vec2(self.x * cos(angle) - self.y * sin(angle), self.x * sin(angle) + self.y * cos(angle))
	
	@staticmethod
	def zero() -> Vec2:
		return Vec2(0, 0)

	@staticmethod
	def distance(a, b) -> float:
		if type(a) != Vec2 or type(b) != Vec2: return 0
		return ((a.x - b.x) ** 2 + (a.y - b.y) ** 2) ** .5

class CosmicBody:
	def __init__(self, mass, radius, airDensity, position):
		self.mass = mass
		self.radius = radius
		self.airDensity = airDensity
		self.position = position

Earth = CosmicBody(5.9726 * 10**24, 6_378_000, getAirDensity, Vec2(0, -6_378_000))
Moon = CosmicBody(7.35 * 10**22, 1_738_000, lambda: 0, Vec2(0, 384_400_000))
Sun = CosmicBody(1.989 * 10**30, 700_000_000, lambda: 0, Vec2(0, 1500_000_000_000))

class RocketStage:
	def __init__(self, dryMass, fuelMass, reserveFuelMass, fuelConsumtion, T):
		self.dryMass = dryMass
		self.fuelMass = fuelMass
		self.reserveFuelMass = reserveFuelMass
		self.fuelConsumtion = fuelConsumtion
		self.T = T
		self.mass = lambda: self.dryMass + self.fuelMass + self.reserveFuelMass

class Rocket:
	def __init__(self, position: Vec2 = Vec2.zero(), speed: Vec2 = Vec2.zero(), diameter: float = 9, cd: float = 0.05):
		self.diameter = diameter
		self.cd = cd
		self.area = self.diameter ** 2 * pi / 4
		self.position = position
		self.thrust = 1
		self.altitude = Vec2.distance(self.position, Earth.position) - Earth.radius
		self.angle = 0
		self.speed = speed
		self.mass = lambda: (sum(map(lambda x: x.mass(), self.stages)) if self.stages else 1)
		
		self.stages = [RocketStage(100_000, 1_200_000, 2_000, 2_000, 3 * 2.2555 * 10 ** 6), RocketStage(200_000, 3_400_000, 50_000, 22_000, 75.9 * 10 ** 6)]

	def update(self, deltaTime):
		#if len(self.stages) <= 0 or len(self.stages) == 1 and self.stages[-1].fuelMass <= 0: system.stop("No active rocket stages left"); return
		if self.altitude < 0: system.stop("Rocket crushed"); return

		self.angle = getAngle(system.time)
		self.thrust = getThrust(system.time)

		EarthDir = (Earth.position - self.position).normalized()
		MoonDir = (Moon.position - self.position).normalized()

		gravityForce = EarthDir * G * self.mass() * Earth.mass / (Vec2.distance(Earth.position, self.position) ** 2)
		dragForce = self.speed.normalized().negative() * Earth.airDensity(self.altitude) * (self.speed ** 2) * self.cd * self.area / 2

		self.speed += ((EarthDir.negative() * self.thrust * (self.stages[-1].T if self.stages and self.stages[-1].fuelMass > 0 else 0)).rotated(-self.angle * pi / 180) + gravityForce + dragForce) * deltaTime / self.mass()
		self.position += self.speed * deltaTime
		self.altitude = Vec2.distance(self.position, Earth.position) - Earth.radius

		if self.stages and self.stages[-1].fuelMass > 0:
			self.stages[-1].fuelMass -= self.stages[-1].fuelConsumtion * deltaTime * self.thrust
		elif self.stages and self.stages[-1].fuelMass < 0:
			self.stages[-1].fuelMass = 0
		if len(self.stages) > 1 and self.stages[-1].fuelMass <= 0: self.stages = self.stages[:-1]
		
	def printReport(self):
		print("Time:", f'{system.time:10.2f}', "|", "Rocket mass:", f'{self.mass(): 10.0f}',\
		"|", "Current rocket stage fuel mass:", f'{self.stages[-1].fuelMass if self.stages else 0:10.0f}',\
		"|", "Speed:", f'{self.speed.x:8.0f}, {self.speed.y:8.0f}', "|", "Altitude:", f'{self.altitude:10.0f}', "|", "Angle:", f'{self.angle:3.1f}')


class System:
	def __init__(self, deltaTime = np.double(.001), simulationTime = 1000, dataCollectionInterval = 1):
		self.rocket = Rocket()
		self.deltaTime = deltaTime
		self.simulationTime = simulationTime
		self.dataCollectionInterval = dataCollectionInterval
		self.maxIterator = self.simulationTime / self.deltaTime
		self.dataCollectionIterations = self.dataCollectionInterval / self.deltaTime
		self.iterator = 0
		self.time = np.double(0)
		self.active = False
	
	def start(self):
		self.active = True
		self.collectData()
		while self.update(): pass
	
	def update(self):
		if not self.active: return False

		self.rocket.update(self.deltaTime)
		self.iterator += 1
		if self.iterator >= self.maxIterator: return False
		if self.iterator % self.dataCollectionIterations == 0:
			self.collectData()

		self.time += self.deltaTime

		return True
	
	def collectData(self):
		global xpoints, ypoints_altitude, ypoints_speedX, ypoints_speedY, ypoints_fuel, xpoints_position, ypoints_position, zpoints_position, ypoints_localSpeedX, ypoints_localSpeedY
		rocket = self.rocket
		dir = rocket.position - Earth.position
		angle = atan2(dir.y, dir.x)
		xpoints = np.append(xpoints, [self.time])
		ypoints_altitude = np.append(ypoints_altitude, [rocket.altitude])
		ypoints_speedX = np.append(ypoints_speedX, [rocket.speed.x])
		ypoints_speedY = np.append(ypoints_speedY, [rocket.speed.y])
		localSpeed = self.rocket.speed.rotated(-angle + pi/2)
		ypoints_localSpeedX = np.append(ypoints_localSpeedX, [localSpeed.x])
		ypoints_localSpeedY = np.append(ypoints_localSpeedY, [localSpeed.y])
		ypoints_fuel = np.append(ypoints_fuel, [rocket.stages[-1].fuelMass / 1_000 if rocket.stages else 0])
		xpoints_position = np.append(xpoints_position, [rocket.position.x])
		ypoints_position = np.append(ypoints_position, [rocket.position.y])
		zpoints_position = np.append(zpoints_position, [self.time])
		self.rocket.printReport()
	
	def stop(self, message = "System stoped"):
		self.collectData()
		self.active = False
		print(f'System has been stoped with message: "{message}"')

deltaTime = np.double(.1)
simulationTime = 10000
dataCollectionInterval = 10

system = System(deltaTime, simulationTime, dataCollectionInterval)

''' # Draw Sun
for i in range(1000):
	pos = Sun.position + Vec2(cos(i / 500 * pi) * Sun.radius, sin(i / 500 * pi) * Sun.radius)
	xpoints_position = np.append(xpoints_position, [pos.x])
	ypoints_position = np.append(ypoints_position, [pos.y])
	zpoints_position = np.append(zpoints_position, [-6000])
'''

# Draw Earth
for i in range(1000):
	pos = Earth.position + Vec2(cos(i / 500 * pi) * Earth.radius, sin(i / 500 * pi) * Earth.radius)
	xpoints_position = np.append(xpoints_position, [pos.x])
	ypoints_position = np.append(ypoints_position, [pos.y])
	zpoints_position = np.append(zpoints_position, [-4000])

''' # Draw Moon
for i in range(1000):
	pos = Moon.position + Vec2(cos(i / 500 * pi) * Moon.radius, sin(i / 500 * pi) * Moon.radius)
	xpoints_position = np.append(xpoints_position, [pos.x])
	ypoints_position = np.append(ypoints_position, [pos.y])
	zpoints_position = np.append(zpoints_position, [-2000])
'''

system.start()

print(f'Max altitude: {np.max(ypoints_altitude):.2f}, min altitude: {np.min(ypoints_altitude[min(int(1000 / dataCollectionInterval), len(ypoints_altitude) - 1):]):.2f}')

plt.figure(1, (12, 8))
plt.subplot(2, 3, 1)
plt.plot(xpoints, ypoints_altitude)
plt.title("Altitude")
plt.subplot(2, 3, 2)
plt.plot(xpoints, ypoints_localSpeedX)
plt.title("Local Speed X")
plt.subplot(2, 3, 3)
plt.plot(xpoints, ypoints_localSpeedY)
plt.title("Local Speed Y")
plt.subplot(2, 3, 4)
plt.plot(xpoints, ypoints_fuel)
plt.title("Fuel")
plt.subplot(2, 3, 5)
plt.plot(xpoints, ypoints_speedX)
plt.title("Speed X")
plt.subplot(2, 3, 6)
plt.plot(xpoints, ypoints_speedY)
plt.title("Speed Y")

plt.figure(2, (10, 8))
plt.scatter(xpoints_position, ypoints_position, c=zpoints_position, cmap='viridis')
plt.colorbar()
plt.title("Scene")

plt.show()

