

#include<iostream>
#include<vector>
#include<cmath>
#include<functional>

using namespace std;

class System;

const double G = 6.6743 * pow(10, -11);

double airDensitySpreadsheet[70] = {
	0, 1.2250, 500, 1.1673, 1000, 1.1117, 1500, 1.0581, 2000, 1.0065, 2500, 0.9569, 3000, 0.9093, 4000, 0.8194, 5000, 0.7365,
	6000, 0.6601, 7000, 0.59, 8000, 0.5258, 9000, 0.4671, 10'000, 0.4135, 11'000, 0.3648, 12'000, 0.3119, 14'000, 0.2279,
	16'000, 0.1665, 18'000, 0.1216, 20'000, 0.0889, 24'000, 0.0469, 28'000, 0.0251, 32'000, 0.0136,
	36'000, 7.26 * pow(10, -3), 40'000, 4.00 * pow(10, -3), 50'000, 1.03 * pow(10, -3),
	60'000, 3.00 * pow(10, -4), 80'000, 1.85 * 10 * pow(10, -5), 100'000, 5.55 * 10 * pow(10, -7),
	150'000, 2.00 * pow(10, -9), 200'000, 2.52 * 10 * pow(10, -10), 300'000, 1.92 * 10 * pow(10, -11),
	500'000, 5.21 * 10 * pow(10, -13), 700'000, 3.07 * pow(10, -14), 1'000'000, 3.56 * pow(10, -15)
};

double getAirDensity(double altitude) {
	int len = sizeof(airDensitySpreadsheet) / sizeof(*airDensitySpreadsheet);
	if (altitude <= airDensitySpreadsheet[0]) return airDensitySpreadsheet[1];
	if (altitude >= airDensitySpreadsheet[len - 2]) return airDensitySpreadsheet[len - 1];
	for(size_t i = 0; i < sizeof(airDensitySpreadsheet); i += 2) {
		if (altitude <= airDensitySpreadsheet[i]) {
			double p[2] = {airDensitySpreadsheet[i-2], airDensitySpreadsheet[i-1]};
			double k = (altitude - p[0]) / (airDensitySpreadsheet[i] - p[0]);
			return p[1] * (1 - k) + airDensitySpreadsheet[i + 1] * k;
		}
	}
	return 0.;
}

// double angleSpreadsheet[24] = {5, 0, 20, 10, 60, 15, 150, 20, 160, 5, 300, 10, 400, 25, 540, 40, 550, 55, 800, 65, 999, 70, 1000, 90};
double angleSpreadsheet[16] = {0, 0, 100'000, 10, 200'000, 20, 300'000, 30, 400'000, 40, 500'000, 50, 600'000, 70, 700'000, 90};

double getAngle(double altitude) {
	int len = sizeof(angleSpreadsheet) / sizeof(*angleSpreadsheet);
	if (altitude <= angleSpreadsheet[0]) return angleSpreadsheet[1];
	if (altitude >= angleSpreadsheet[len - 2]) return angleSpreadsheet[len - 1];
	for(size_t i = 0; i < sizeof(angleSpreadsheet); i += 2) {
		if (altitude <= angleSpreadsheet[i]) {
			double p[2] = {angleSpreadsheet[i-2], angleSpreadsheet[i-1]};
			double k = (altitude - p[0]) / (angleSpreadsheet[i] - p[0]);
			return p[1] * (1 - k) + angleSpreadsheet[i + 1] * k;
		}
	}
	return 0;
}

//double throttlingSpreadsheet[10] = {0, 1, 600, 1, 601, .2, 2000, .2, 2001, 0};
double throttlingSpreadsheet[10] = {0, 1, 999999999, 1, 601, .2, 2000, .2, 2001, 0};

double getThrottling(double time) {
	if (time <= throttlingSpreadsheet[0]) return throttlingSpreadsheet[1];
	if (time >= throttlingSpreadsheet[8]) return throttlingSpreadsheet[9];
	for(size_t i = 0; i < sizeof(throttlingSpreadsheet); i += 2) {
		if (time <= throttlingSpreadsheet[i]) {
			double p[2] = {throttlingSpreadsheet[i-2], throttlingSpreadsheet[i-1]};
			double k = (time - p[0]) / (throttlingSpreadsheet[i] - p[0]);
			return p[1] * (1 - k) + throttlingSpreadsheet[i + 1] * k;
		}
	}
	return 0;
}

class Vec2 {
	public:
		double x;
		double y;

		Vec2(double x, double y) : x(x), y(y) {}

		double length() {
			return pow(pow(x, 2) + pow(y, 2), .5);
		}

		Vec2 normalized() {
			return *this / length();
		}

		Vec2 negative() {
			return *this * (-1);
		}

		Vec2 rotated(double angle) {
			return Vec2(x * cos(angle) - y * sin(angle), x * sin(angle) + y * cos(angle));
		}
	
	Vec2 operator + (const Vec2& other) const {
		return Vec2(x + other.x, y + other.y);
	}
	
	Vec2 operator + (const double& other) const {
		return Vec2(x + other, y + other);
	}
	
	Vec2& operator += (const Vec2& other) {
		x += other.x;
		y += other.y;
		return *this;
	}
	
	Vec2 operator - (const Vec2& other) const {
		return Vec2(x - other.x, y - other.y);
	}
	
	Vec2 operator - (const double& other) const {
		return Vec2(x - other, y - other);
	}
	
	Vec2& operator -= (const Vec2& other) {
		x -= other.x;
		y -= other.y;
		return *this;
	}
	
	Vec2 operator * (const Vec2& other) const {
		return Vec2(x * other.x, y * other.y);
	}
	
	Vec2 operator * (const double& other) const {
		return Vec2(x * other, y * other);
	}
	
	Vec2& operator *= (const Vec2& other) {
		x *= other.x;
		y *= other.y;
		return *this;
	}
	
	Vec2 operator / (const Vec2& other) const {
		if (other.x == 0 || other.y == 0) return Vec2::zero();
		return Vec2(x / other.x, y / other.y);
	}
	
	Vec2 operator / (const double& other) const {
		if (other == 0.) return Vec2::zero();
		return Vec2(x / other, y / other);
	}
	
	Vec2& operator /= (const Vec2& other) {
		if (other.x == 0 || other.y == 0) {
			x = y = 0;
			return *this;
		}
		x /= other.x;
		y /= other.y;
		return *this;
	}

	static Vec2 zero() {
		return Vec2(0., 0.);
	}

	static double distance(Vec2 a, Vec2 b) {
		return (a - b).length();
	}
};

struct CosmicBody {
	double mass;
	double radius;
	Vec2 position;

	function<double(double altitude)> getAirDensity;
};

CosmicBody Earth{5.9726 * pow(10, 24), 6'378'000, Vec2(0, -6'378'000), getAirDensity};
CosmicBody Moon{7.35 * pow(10, 22), 1'738'000, Vec2(0, -384'400'000)};
CosmicBody Sun{1.989 * pow(10, 30), 700'000'000, Vec2(0, 1'500'0000'000'000)};

class RocketEngine {
	private:
		double thrustSeaLevel;
		double thrustVacuum;
		double fuelConsumptionSeaLevel;
		double fuelConsumptionVacuum;

	public:
		double getThrust(double airDensity) {
			double k = min(airDensity, 1.);
			return thrustSeaLevel * k + thrustVacuum * (1 - k);
		}

		double getFuelConsumption(double airDensity) {
			double k = min(airDensity, 1.);
			return fuelConsumptionSeaLevel * k + fuelConsumptionVacuum * (1 - k);
		}

		RocketEngine(double fuelConsumptionSeaLevel, double fuelConsumptionVacuum, double thrustSeaLevel, double thrustVacuum) :
		fuelConsumptionSeaLevel(fuelConsumptionSeaLevel), fuelConsumptionVacuum(fuelConsumptionVacuum), thrustSeaLevel(thrustSeaLevel), thrustVacuum(thrustVacuum) {}
};

class RocketStage {
	private:
		RocketEngine* rocketEngine;
		int enginesQuantity;

	public:
		double dryMass;
		double fuelMass;
		double reserveFuelMass;

		RocketStage(double dryMass, double fuelMass, double reserveFuelMass, RocketEngine* rocketEngine, int enginesQuantity) :
		dryMass(dryMass), fuelMass(fuelMass), reserveFuelMass(reserveFuelMass), rocketEngine(rocketEngine), enginesQuantity(enginesQuantity) {}

		double getMass() {
			return dryMass + fuelMass + reserveFuelMass;
		}

		double getThrust(double airDensity) {
			return rocketEngine->getThrust(airDensity) * enginesQuantity;
		}
		
		double getFuelConsumption(double airDensity) {
			return rocketEngine->getFuelConsumption(airDensity) * enginesQuantity;
		}
};

RocketEngine* Merlin1D = new RocketEngine(300, 280, 845'000, 914'000);
RocketEngine* Merlin1DVac = new RocketEngine(300, 280, 845'000, 981'000); // Specific impulse in vacuum = 348 s

class Rocket {
	private:
		System* system;
		vector<RocketStage*> stages;
		double altitude = 1000;
		Vec2 position = Vec2::zero();
		Vec2 speed = Vec2::zero();
		double diameter = 9;
		double cd = .5;
		double throttling = 1;
		double angle = 0;
		double area;
		double firstCosmicVelocity = 0;

	public:
		Rocket(System* system) : system(system) {
			stages.push_back(new RocketStage(4'000, 110'000, 0, Merlin1DVac, 1));
			stages.push_back(new RocketStage(22'000, 430'000, 0, Merlin1D, 9));
			area = pow(diameter, 2);
			position += Vec2(0, altitude);
		}

		void update(float deltaTime);

		double getMass() {
			double mass = 0;
			for(size_t i = 0; i < stages.size(); i++) {
				mass += stages[i]->getMass();
			}
			mass += 100;
			return mass;
		}

		void printReport();
};

class System {
	private:
		bool active = false;
		int iterator = 0;
		Rocket* rocket;
		int maxIterator;
		int dataCollectionIterations;
		double dataCollectionInterval;
		double simulationTime;
		double deltaTime;
		
		bool update() {
			if (!active) return false;
			if (iterator++ > maxIterator) { stop("System completed successfully"); return false; }
			if (iterator % dataCollectionIterations == 0) collectData();
			time += deltaTime;
			rocket->update(deltaTime);
			return true;
		}

		void collectData() {
			rocket->printReport();
		}
	
	public:
		double time = 0;

		System(double deltaTime = .1, double simulationTime = 1000, double dataCollectionInterval = 10) :
		deltaTime(deltaTime), simulationTime(simulationTime), dataCollectionInterval(dataCollectionInterval) {
			rocket = new Rocket(this);
			maxIterator = simulationTime / deltaTime;
			dataCollectionIterations = dataCollectionInterval / deltaTime;
		}

		void start() {
			collectData();
			active = true;
			while(update()) {}
		}

		void stop(string message = "...") {
			collectData();
			active = false;
			printf("System has been stoped with message: \"%s\"\n", message.c_str());
		}
};

void Rocket::update(float deltaTime) {
	if (altitude < 0.) system->stop("Rocket has been crushed");

	angle = getAngle(altitude);
	throttling = getThrottling(system->time);

	double airDensity = Earth.getAirDensity(altitude);

	Vec2 EarthDir = (Earth.position - position).normalized();
	Vec2 MoonDir = (Moon.position - position).normalized();

	Vec2 gravityForce = EarthDir * G * getMass() * Earth.mass / pow(Vec2::distance(Earth.position, position), 2);

	Vec2 dragForce = speed.normalized().negative() * airDensity * pow(speed.length(), 2) * cd * area / 2;
	
	firstCosmicVelocity = sqrt(2 * G * Earth.mass / (Earth.radius + altitude));
	
	speed += ((EarthDir.negative() * throttling * (stages.empty() || stages.back()->fuelMass <= 0 ? 0 : stages.back()->getThrust(airDensity))).rotated(-angle * 3.14159265359 / 180) + gravityForce + dragForce) * deltaTime / getMass();

	position += speed * deltaTime;

	altitude = Vec2::distance(position, Earth.position) - Earth.radius;

	if (stages.back()->fuelMass > 0) stages.back()->fuelMass -= stages.back()->getFuelConsumption(airDensity) * deltaTime * throttling;
	else stages.back()->fuelMass = 0;

	if (stages.size() > 1 && stages.back()->fuelMass <= 0) stages.pop_back();
}

void Rocket::printReport() {
	printf("Time: %7.2f  |  Rocket mass: %9.0f  |  Current rocket stage fuel mass: %9.0f  |  Speed: (%12.2f, %12.2f)  |  Position: (%14.2f, %14.2f)  |  Angle: %03.1f  |  Altitude: %11.0f  |  %10f  %10f\n", system->time, getMass(), stages.back()->fuelMass, speed.x, speed.y, position.x, position.y, angle, altitude, speed.length(), firstCosmicVelocity);
}

int main() {
	double deltaTime = .1;
	double simulationTime = 10000;
	double dataCollectionInterval = 10;

	int setAngles = 0;

	int x = 1;
	while (x)
	{
		printf("Enter params (dt, t, interval, setangles):\n");
		scanf("%lf %lf %lf %i", &deltaTime, &simulationTime, &dataCollectionInterval, &setAngles);
		
		if (setAngles == 1) {
			angleSpreadsheet[0] = 0;
			angleSpreadsheet[2] = 2'000;
			angleSpreadsheet[4] = 10'000;
			angleSpreadsheet[6] = 50'000;
			angleSpreadsheet[8] = 150'000;
			angleSpreadsheet[10] = 300'000;
			angleSpreadsheet[12] = 500'000;
			angleSpreadsheet[14] = 700'000;
			scanf("%lf %lf %lf %lf %lf %lf %lf %lf", &angleSpreadsheet[1], &angleSpreadsheet[3], &angleSpreadsheet[5],
				&angleSpreadsheet[7], &angleSpreadsheet[9], &angleSpreadsheet[11], &angleSpreadsheet[13], &angleSpreadsheet[15]);
		}

		System* system = new System(deltaTime, simulationTime, dataCollectionInterval);

		system->start();

		scanf("%i", &x);
	}

	return 0;

}