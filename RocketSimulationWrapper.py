
import subprocess
import time
import matplotlib.pyplot as plt
import numpy as np
from math import pi, cos, sin
from scipy.optimize import minimize

startTime = time.time()

def main(deltaTime = 1, simulationTime = 1e4, dataCollectionInterval = 10, showGraph: bool = False, printOutput: bool = False, maxDataTime: float = 1e10, angles = [0, 10, 25, 45, 70, 85, 90, 100]):
    result = subprocess.run(["/Programming/RocketSimulation.exe"], input=f'{deltaTime} {simulationTime} {dataCollectionInterval} 1\n{" ".join(map(str, angles))} \n0', capture_output=True, text=True)

    sts = result.stdout

    h = sts.split("\n")

    minMaxAltitude = [10e10, 0]

    if printOutput:
        h[0] = h[0][14:]
        print(*h[:min(len(h), int(maxDataTime / dataCollectionInterval))], sep="\n")

    sts = sts.replace("\n", " ").replace("  ", " ").replace("  ", " ").replace("  ", " ").replace("  ", " ")

    tt = sts.split(" ")

    g = " ".join(h[min(len(h)-1, int(1_000 / dataCollectionInterval)):]).replace("\n", " ").replace("  ", " ").replace("  ", " ").replace("  ", " ").replace("  ", " ").split(" ")

    for i, el in enumerate(g):
        if el == "Altitude:":
            altitude = float(g[i + 1])
            minMaxAltitude[0] = min(minMaxAltitude[0], altitude)
            minMaxAltitude[1] = max(minMaxAltitude[1], altitude)

    if showGraph:
        t = []
        xs = np.array([])
        ys = np.array([])
        times = np.array([])

        for i, el in enumerate(tt):
            if el == "Time:":
                times = np.append(times, np.double(tt[i + 1]))
            if el == "Position:":
                xs = np.append(xs, np.double(tt[i + 2][:-1] if tt[i + 2][-1] == ',' else tt[i + 2]))
                ys = np.append(ys, np.double(tt[i + 3][:-1] if tt[i + 3][-1] == ')' else tt[i + 3]))

        for i in range(1000):
            pos = [cos(i / 500 * pi) * 6_378_000, sin(i / 500 * pi) * 6_378_000 - 6_378_000]
            xs = np.append(xs, [np.double(pos[0])])
            ys = np.append(ys, [np.double(pos[1])])
            times = np.append(times, [-simulationTime/2])

        plt.figure(2, (10, 8))
        plt.scatter(xs, ys, c=times, cmap='viridis')
        plt.colorbar()
        plt.title("Scene")

        plt.show()

    if printOutput:
        print(f'Min altitude: {minMaxAltitude[0]}')
        print(f'Max altitude: {minMaxAltitude[1]}')
        print(f'Time: {time.time() - startTime}')
    
    success = h[-2] != "System has been stoped with message: \"Rocket has been crushed\""

    return {"output": sts, "minMaxAltitude": minMaxAltitude, "success": success}

def f(out):
    res = 0
    if not out["success"]: res += 1e20
    res += ((out["minMaxAltitude"][1] - out["minMaxAltitude"][0]) / 1e2) ** 2
    res += (1e6 / out["minMaxAltitude"][1]) ** 2
    return res

def ff(angles):
    out = main(angles=angles)
    return f(out)

def fMaxOrb(out, orb = 1e5):
    res = 0
    if not out["success"]: res += 1e20
    res += ((out["minMaxAltitude"][1] - orb) / 10) ** 2
    return res

def ffMaxOrb(angles):
    out = main(angles=angles)
    return fMaxOrb(out)

def fOrbDiff(out):
    res = 0
    if not out["success"]: res += 1e20
    res += ((out["minMaxAltitude"][1] - out["minMaxAltitude"][0]) / 1e2) ** 2
    return res

def ffOrbDiff(angles):
    out = main(angles=angles)
    return fOrbDiff(out)

def ffnew(angles):
    out = main(angles=angles)
    return fOrbDiff(out)

# --- Simulate flight ---
# main(printOutput=True, showGraph=True, maxDataTime=10, angles=[0.0, 0.07, 3.73, 148.09, 102.39, 92.84, 80.42, 57.0])
# main(printOutput=True, showGraph=True, dataCollectionInterval=10, maxDataTime=1000, angles=[0, 50, 90, 100, 110, 110, 90, 90])
# main(printOutput=True, showGraph=True, dataCollectionInterval=10, maxDataTime=1000, angles=[-0.0, 2.21, 3.95, 42.24, 74.63, 184.34, -28.06, 171.66])
# main(printOutput=True, showGraph=True, dataCollectionInterval=10, maxDataTime=1000, angles=[-0.0, -0.76, 2.15, 44.67, 67.86, 178.71, -95.39, 213.93])
# main(printOutput=True, showGraph=True, dataCollectionInterval=10, maxDataTime=1000, angles=[0.05, -41.54, 3.71, 117.19, 47.51, 152.01, -316.94, -1267.4])
# main(printOutput=True, showGraph=True, dataCollectionInterval=10, maxDataTime=1000, angles=[0.0, -0.76, 2.19, 44.71, 67.88, 178.7, -95.23, 216.33])
main(printOutput=True, showGraph=True, dataCollectionInterval=10, maxDataTime=1000, angles=[0.0, -0.76, 2.2, 44.73, 67.87, 178.7, -96.69, 218.92])


# --- Check angles' error ---
# out = main(angles=[0.0, 8.17, 19.95, 39.72, 45.87, 67.37, 109.33, 137.88])
# print(f(out))


# --- Angles finding method ---
def findAnglesMethod(angles, fun=ff, iters=100):
    return minimize(fun, np.array(angles), method="Nelder-Mead", options={"maxiter": iters, "disp": True})


# --- Test angles finding method ---
def testFindAnglesMethod(angles):
    #angles = [0.0, 8.17, 19.95, 39.72, 45.87, 67.37, 109.33, 137.88]
    for _ in range(3):
        res = minimize(ffMaxOrb, np.array(angles), method="Nelder-Mead", options={"maxiter": 100, "disp": True})
        angles = list(map(lambda x: float(f'{float(x):5.2f}') if x != '' else '', str(res.x).replace("[", "").replace("]", "").replace("  ", " ").split(" ")))
        print(angles)
        res = minimize(ffOrbDiff, np.array(angles), method="Nelder-Mead", options={"maxiter": 100, "disp": True})
        angles = list(map(lambda x: float(f'{float(x):5.2f}') if x != '' else '', str(res.x).replace("[", "").replace("]", "").replace("  ", " ").split(" ")))
        print(angles)
    return res


def findAngles():
    #res = testFindAnglesMethod([0.0, 50.0, 90.0, 100.0, 110.0, 110.0, 90.0, 90.0])
    res = findAnglesMethod([0.0, -0.76, 2.19, 44.71, 67.88, 178.7, -95.23, 216.33], ff, 200)
    angles = list(map(lambda x: float(f'{float(x):5.2f}') if x != '' else '', (str(res.x)[2:] if str(res.x)[1] == " " else str(res.x)).replace("[", "").replace("]", "").replace("  ", " ").split(" ")))
    print("Raw angles:", res.x)
    print("Angles:", angles)
    print("Funciton value:", res.fun)
    main(maxDataTime=10, angles=angles, printOutput=True)

# --- No comments ---
# findAngles()

print("Time:", time.time() - startTime)

