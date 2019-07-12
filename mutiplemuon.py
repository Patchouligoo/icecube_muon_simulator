from math import *
import numpy as np
import sys
import os
import subprocess
import time
from scipy.interpolate import UnivariateSpline


def convertcord(x, y, z):
    R = sqrt(pow(x, 2) + pow(y,2) + pow(z,2))
    theta = acos(z/R)
    try:
        phi = atan(y/x)
    except ZeroDivisionError:
        phi = pi/2
    if(x < 0 and y > 0):
        phi = pi + phi
    if(x <= 0 and y < 0):
        phi = phi - pi

    depth = R/1000 - 6371

    latitude = degrees(theta) - 90  # since we use south pole as theta = 0, north pole as theta = 180

    longitude = -1 * degrees(phi)

    return depth, latitude, longitude


def getB(x, y, z):

    depth1, latitude1, longitude1 = convertcord(x, y, z)

    main = "./wmm_point.exe"

    latitude = latitude1
    longitude = longitude1
    height = depth1
    date = 2019.5

    p = subprocess.Popen(main, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)


    if(height < -10):
        rs, out = p.communicate(os.linesep.join(["c", str(latitude), str(longitude), str(height), "c", str(date), "c", "c"]))
    else:
        rs, out = p.communicate(os.linesep.join(["c", str(latitude), str(longitude), str(height), str(date), "c", "c"]))


    #print(str(rs))


    pre = str(rs).split('Secular Change')[1]

    array = pre.split('+/-')

    array = [x for x in array if x != ' ']

    #print(array)

    B_theta = float(array[2].split('X\t=\t')[1])
    B_phi = -1 * float(array[3].split('Y\t=\t')[1])
    B_r = -1 * float(array[4].split('Z\t=\t')[1])

    R = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))
    try:
        theta = acos(z / R)
    except ZeroDivisionError:
        theta = 0

    try:
        phi = atan(y / x)
    except ZeroDivisionError:
        phi = pi/2
    if (x < 0 and y > 0):
        phi = pi + phi
    if (x <= 0 and y < 0):
        phi = phi + pi
    if (x > 0 and y < 0):
        phi = phi + 2 * pi


    thetap = theta + pi / 2
    phip = phi + pi / 2

    B_x = B_theta * sin(thetap) * cos(phi) + B_phi * cos(phip) + B_r * sin(theta) * cos(phi)

    B_y = B_theta * sin(thetap) * sin(phi) + B_phi * sin(phip) + B_r * sin(theta) * sin(phi)

    B_z = B_r * cos(theta) + B_theta * cos(thetap)

    return np.array([B_x, B_y, B_z]) * pow(10, -9)



class muon(object):
    v_dir = np.array([0, 0, 0])
    """ Note that input v_dir is the direction we received from ice cube detector """
    def __init__(self, v_d, initial_E):
        self.initialvdir = v_d/np.linalg.norm(v_d) * -1
        self.v_dir = v_d * -1
        self.v_dir = self.v_dir / np.linalg.norm(self.v_dir)
        self.initialE = initial_E
        self.E = initial_E



    timedelta1 = 0.0

    timedelta2 = 0.0


    updateref = 5 * 10e-8

    finalvdir = np.array([0,0,0])


    const_c = 299792458  # m/s
    const_mass = 105.658 * np.power(10, 6)  # ev/c^2
    const_e = 1.602 * pow(10, -19)

    Radius = 0



    """
    Get E loss =======================================================================
    """

    """
    brems loss
    """
    x = np.log10(np.array([1.461684024050507, 2.5826583413126927, 4.53013361076875, 6.767498878893491, 14.938550176741655, 24.898127176298612, 45.13080182813426, 78.58645531599923, 145.06998863954541, 258.20249675073455, 494.3531785103597, 919.2562918101086, 1642.1172873211872, 2998.3335552024096, 5740.59407597251, 10597.09226019805, 20289.138593628002, 38004.261024316555, 72762.76360892662, 132857.1578792675, 251598.32534574467, 476464.4850698085, 902304.9148712363, 1665647.8934370114, 3154321.7339125252, 5973498.744985013, 11189159.386201639, 21422707.28221403, 41015805.69719786, 76828071.26806813, 147094631.95801523, 268579336.0773514, 508622284.70174766, 963203767.920256, 1824067734.4285278, 3367214928.4157166, 6541638461.111238, 12075813999.008741, 22868580269.4823, 41755623737.25711, 68585830152.334236, 161730231721.6129, 372759372031.49225, 926900181144.2858, 2136339791152.6511, 6670767806717.104, 22472342828110.438, 75704357701613.97, 219108557989221.94, 634158477023743.6]))

    y = np.log10(np.array([3.733156987116539e-07, 8.095239915988781e-07, 1.7078461796892677e-06, 2.9052353010850296e-06, 7.786269768517139e-06, 1.490323867310003e-05, 3.0310080833483925e-05, 5.883981279000656e-05, 0.000119759295481561, 0.00023480269647880865, 0.00047030464919945846, 0.0009448905664426858, 0.0017767448378974262, 0.003406608103449301, 0.006623223148092696, 0.012575918623439807, 0.024864526959477308, 0.04656754157363181, 0.08964392495469672, 0.16721819584780295, 0.313713021580771, 0.5988574584781586, 1.1493041266144675, 2.136106751642162, 3.991457647717873, 7.7012924336549045, 14.13993472526518, 27.17827773683473, 52.63947912275994, 99.5690475321118, 190.65197959335754, 349.97919737975246, 656.5854328457507, 1251.705076910184, 2367.18339199842, 4423.239731010238, 8500.258391888041, 15947.078790179039, 30158.58987820563, 55298.6685535126, 94500.9791865196, 263444.81485700083, 621016.9418915603, 1610262.027560936, 4592709.387993489, 11908637.320879903, 33965223.9733046, 117210229.75334746, 367719991.4775928, 1268961003.1679182]))

    spl_brems = UnivariateSpline(x, y)

    #print(spl_brems(E))


    """
    ionization
    """

    x2 = np.log10(np.array([1.5220674559482037, 2.8939837612399706, 5.5024775528776875, 10.462138601271949, 19.892192755058282, 37.82203120079919, 71.91293899916963, 136.73170454655556, 259.9748986538075, 494.30341085995195, 939.844051302618, 1786.9729833185027, 3397.662025614108, 6460.146486860787, 12283.002934689117, 23354.294117082045, 44404.69945397581, 84428.89875895283, 160528.9312460653, 305221.768207315, 580333.5701824363, 1103417.540166852, 2097983.5227610967, 3989002.079042055, 7584491.209759514, 14420776.367390024, 27418950.762403015, 52133036.513286345, 99123176.50835826, 188467904.00561142, 358343548.819518, 681336695.8054671, 1295459886.4145386, 2463123339.225409, 4683260861.923352, 8904520513.260454, 16930614780.768806, 32191033355.240383, 61206438271.53116, 86898613578.95168]))

    y2 = np.log10(np.array([0.002278580087571453, 0.002418516827312156, 0.0025332715618363102, 0.0026534712239900744, 0.0027428043864998015, 0.0028164315316423286, 0.0028920351051909504, 0.0029795176895412437, 0.0030900425993135735, 0.0031312421028716012, 0.003204667419525443, 0.0032581658822394195, 0.0033678570142216483, 0.003424079751034714, 0.0035628735997951092, 0.003598442364569703, 0.0037195894241758612, 0.0037942267022630866, 0.003870361651906195, 0.003974256658236785, 0.004121681310255326, 0.004162828746136487, 0.004317248460837103, 0.004389320262859955, 0.004537093435149979, 0.004643484928116422, 0.004783947974769667, 0.004879942670589285, 0.00506096365437274, 0.005179639517396888, 0.005318680446604392, 0.005407470074453622, 0.005571043281452727, 0.005720590954009233, 0.0059327957015278005, 0.006011897600533528, 0.006214297317735928, 0.006402276636624847, 0.006487638119708833, 0.00668388585820374]))

    spl_ionization = UnivariateSpline(x2, y2)

    """
    photo
    """
    x3 = np.log10(np.array([1.222371229325602, 2.0477641356023994, 4.112599491740337, 7.049800885528415, 14.438361174398503, 25.54227594473224, 43.672769627443735, 82.10396294843815, 165.67642589467212, 296.60543180286766, 561.6967209088427, 1033.111185800716, 1907.115257773746, 3515.385337410887, 6864.473286104272, 12564.33360000109, 23649.451793280514, 43024.13507255546, 82977.13679723056, 160499.0636303774, 285872.45005424676, 525029.3699770753, 992100.2654593317, 1791749.2562108133, 3493635.8422306427, 6449224.741148545, 12302679.298717404, 22777042.983337004, 42695584.67702104, 81685113.76956183, 153903013.77496335, 283689407.7599781, 508622284.70174766, 958295812.7417725, 1778066705.107017, 3367214928.4157166, 6170655653.618712, 11771274938.663185, 22985702733.616886, 41150599717.373604, 71264586069.86523, 128792362543.36287, 275144248897.9444, 684171273157.8615, 1576892900167.8533, 3634457212917.827, 10519086393080.568, 32846077743642.887, 95065290143238.84, 320253957723048.3, 926900181144285.8]))

    y3 = np.log10(np.array([4.4292991010359377e-07, 9.459283316511261e-07, 2.058275523935026e-06, 3.5669072658896587e-06, 6.978425262098202e-06, 1.17808023191723e-05, 1.9106797900419316e-05, 3.493400532934057e-05, 6.988663893751543e-05, 0.00012432353065317033, 0.00023669220700157982, 0.0004437969745912459, 0.0008388138535458333, 0.0015767395566546791, 0.0032023565085742044, 0.005920209554888535, 0.011732001214515765, 0.022276267802683693, 0.04509826086588131, 0.0878554967503327, 0.16390089895919477, 0.31071018587669436, 0.6052915833443382, 1.1440523361195993, 2.253350831083066, 4.347713351214447, 8.50472379202314, 16.447005580681505, 30.979634646219335, 62.018381038931594, 120.04521836528745, 226.89562566189807, 405.91817273497117, 790.7653648036409, 1547.908591959046, 2931.2667980579117, 5574.279032604014, 10681.621457295676, 21764.026724240746, 38976.17026610026, 69132.46177288817, 135217.99213228124, 318748.5925972282, 751384.2180360104, 1948297.047624233, 4592709.387993489, 11908637.320879903, 37360595.771877415, 106558035.900799, 404479533.8205424, 1268961003.1679182]))

    spl_photo = UnivariateSpline(x3, y3)

    """=============================================================================================================="""


    v_total = 0
    pos = np.array([0.0001, 0.0001, 6371000 - 2450])

    prepos = np.array([])
    preenergy = 0
    totalDis = 0
    predis = 0
    totalAng = 0

    def getEloss(self, E):
        if E > pow(10, 28):
            print(E)
            print("!!!!!!!!!!!!!!!!! out of range")
            exit()
        E = E / pow(10, 9)

        totalloss = (pow(10, self.spl_brems(log10(E))) * 2 + pow(10, self.spl_ionization(log10(E))) + pow(10, self.spl_photo(log10(E)))) * 0.9167 * pow(10, 9)

        return totalloss

    def getnewE(self, E, dis):

        deltax = dis/10
        num = 10
        if dis > 1000:
            deltax = dis/100
            num = 100
        if dis > 10000:
            deltax = dis / 500
            num = 500

        for i in range(0, num):
            E_loss1 = self.getEloss(E)
            E_temp = E + 100 * deltax * E_loss1
            E_loss2 = self.getEloss(E_temp)
            E_loss = (E_loss1 + E_loss2)/2

            E = E + 100 * deltax * E_loss

        return E


    def convertev(self, E) :
        return E * self.const_e


    def calV(self, E) :
        temp = self.convertev(self.const_mass) * self.const_c
        E = self.convertev(E)
        temp2 = temp/E
        temp3 = np.power(temp2, 2)
        temp4 = np.power(self.const_c, 2)
        if(temp4 <= temp3): return -1
        return np.power(temp4 - temp3, 0.5)


    def calR(self, v_dir, E, B, cosT):
        temp1 = np.power(self.convertev(self.const_mass), 2)
        temp2 = sqrt(np.power(self.convertev(E), 2) - temp1) * pow(cosT, 2)
        temp3 = self.const_e * self.const_c * np.linalg.norm(np.cross(v_dir, B))
        if temp3 == 0:
            return 0
        return temp2/temp3


    def updatepos(self, t, p):



        """"1. calculate first derivative approximation"""

        const_B = getB(self.pos[0], self.pos[1], self.pos[2]) * -1

        self.v_total = self.calV(self.E)
        if(self.v_total == -1): return


        temp1 = np.cross(const_B, self.v_dir)

        temp2 = np.cross(temp1, const_B)
        # print(temp2)

        if(np.linalg.norm(temp2) == 0):
            v_perpendiculardir = np.array([0,0,0])
            cosT = 0
        else:
            temp3 = temp2/np.linalg.norm(temp2)
            cosT = np.dot(temp3, self.v_dir)
            v_perpendiculardir = cosT * temp3


        v_perpendicular = self.v_total * cosT  # scalar

        v_paralleldir = self.v_dir - v_perpendiculardir

        v_omg = 0

        R = self.calR(self.v_dir, self.E, const_B, cosT)
        if(R == 0):
            v_omg = 0
            omg = 0
        else:
            v_omg = v_perpendicular/R
            omg = v_omg * t



        rotaxis = const_B/np.linalg.norm(const_B)

        M = np.array([[cos(omg) + (1 - cos(omg)) * np.power(rotaxis[0], 2),
                       (1 - cos(omg)) * rotaxis[0] * rotaxis[1] - sin(omg) * rotaxis[2],
                       (1 - cos(omg)) * rotaxis[0] * rotaxis[2] + sin(omg) * rotaxis[1]],
                      [(1 - cos(omg)) * rotaxis[0] * rotaxis[1] + sin(omg) * rotaxis[2],
                       cos(omg) + (1 - cos(omg)) * np.power(rotaxis[1], 2),
                       (1 - cos(omg)) * rotaxis[1] * rotaxis[2] - sin(omg) * rotaxis[0]],
                      [(1 - cos(omg)) * rotaxis[0] * rotaxis[2] - sin(omg) * rotaxis[1],
                       (1 - cos(omg)) * rotaxis[1] * rotaxis[2] + sin(omg) * rotaxis[0],
                       cos(omg) + (1 - cos(omg)) * np.power(rotaxis[2], 2)]])

        old_perp = v_perpendiculardir
        v_perpendiculardir = np.dot(M, v_perpendiculardir)

        v_dir_new = v_perpendiculardir + v_paralleldir

        Dis = self.v_total * t

        E_new = self.getnewE(self.E, Dis)

        v_total_new = self.calV(E_new)

        pos_new = self.pos + v_dir_new * v_total_new * t

        #todo: maintain this? or lose some accuracy??????????????????
        #const_B_new = const_B
        const_B_new = getB(pos_new[0], pos_new[1], pos_new[2]) * -1

        """"2. calculate omg on next[][ step"""

        temp_a = np.cross(v_dir_new, const_B_new)

        temp_b = np.cross(const_B_new, temp_a)

        if (np.linalg.norm(temp_b) == 0):
            # v_perp_new_dir = np.array([0, 0, 0])
            cosT_new = 0
        else:
            temp_c = temp_b / np.linalg.norm(temp_b)
            cosT_new = np.dot(temp_c, v_dir_new)
            # v_perp_new_dir = cosT_new * temp_c

        v_perp_new = v_total_new * cosT_new

        v_omg2 = 0

        R2 = self.calR(v_dir_new, E_new, const_B_new, cosT_new)
        if (R2 == 0):
            v_omg2 = 0
        else:
            v_omg2 = v_perp_new / R2

        """3. calculate averaged v_dir_next as new self.v_dir"""

        v_omg_final = (v_omg + v_omg2) / 2

        omg = v_omg_final * t

        M = np.array([[cos(omg) + (1 - cos(omg)) * np.power(rotaxis[0], 2),
                       (1 - cos(omg)) * rotaxis[0] * rotaxis[1] - sin(omg) * rotaxis[2],
                       (1 - cos(omg)) * rotaxis[0] * rotaxis[2] + sin(omg) * rotaxis[1]],
                      [(1 - cos(omg)) * rotaxis[0] * rotaxis[1] + sin(omg) * rotaxis[2],
                       cos(omg) + (1 - cos(omg)) * np.power(rotaxis[1], 2),
                       (1 - cos(omg)) * rotaxis[1] * rotaxis[2] - sin(omg) * rotaxis[0]],
                      [(1 - cos(omg)) * rotaxis[0] * rotaxis[2] - sin(omg) * rotaxis[1],
                       (1 - cos(omg)) * rotaxis[1] * rotaxis[2] + sin(omg) * rotaxis[0],
                       cos(omg) + (1 - cos(omg)) * np.power(rotaxis[2], 2)]])

        v_perp = np.dot(M, old_perp)

        v_old_dir = self.v_dir

        self.v_dir = v_perp + v_paralleldir

        """4. calculate new position"""

        tempv = (v_old_dir * self.v_total + self.v_dir * v_total_new)

        new_v = tempv / np.linalg.norm(tempv)

        avg_v = self.calV((self.E + E_new) / 2)

        # self.v_dir = new_v

        self.prepos = self.pos

        self.pos = self.pos + new_v * avg_v * t

        Dis = avg_v * t

        self.predis = self.totalDis

        self.totalDis = self.totalDis + Dis

        self.v_total = self.calV(E_new)

        self.totalAng = self.totalAng + omg

        self.finalvdir = self.v_dir

        self.Radius = R

        #E_tempnew_2 = self.E + self.getEloss(self.E) * Dis * 100

        self.preenergy = self.E

        self.E = self.getnewE(self.E, Dis)

        """5. get new update factor"""

        if v_perp_new != 0:
            timet = self.Radius * self.updateref / v_perp_new
        else:
            timet = 10e-5



        #p.write(str(self.pos[0]) + ' ')
        #p.write(str(self.pos[1]) + ' ')
        #p.write(str(self.pos[2]) + '\n')
        #p.write(str(self.omg) + '\n')
        """
        p.write(str(self.totalDis/1000) + ' ')
        p.write(str((6371000 - np.linalg.norm(self.pos))/1000) + '\n')
        """
        #print(self.totalDis)

        if Dis >= 1000:
            timet = 1000/self.v_total

        return timet




def main(v_dir, p, E):
    x = muon(v_dir, E)



    j = 0

    t = 1./10000000000000000


    while 1:
        j = j + 1

        t = x.updatepos(t, p)

        if np.linalg.norm(x.pos) >= 6371000 : #6371000:
            break
        #print(x.pos)

    try:
        defang = acos(np.dot(x.initialvdir, x.finalvdir)/(np.linalg.norm(x.initialvdir) * np.linalg.norm(x.finalvdir)))
    except:
        defang = x.totalAng

    return x.totalAng, defang, x.totalDis, x.E, x.pos, x.prepos, x.predis, x.preenergy



def changecord(theta, phi):
    x = sin(theta) * cos(phi)
    y = sin(theta) * sin(phi)
    z = cos(theta)

    return np.array([x, y, z])




"""
This is the start case for mutiple muons with different initial v_dir
theta change and phi are the same
A plot will be given for each muon in sale figure
"""
def start():


    E = np.array([pow(10, 1), pow(10, 2), pow(10, 3)])* pow(10, 9) #, pow(10, 1.5), pow(10, 2), pow(10, 2.5), pow(10, 3), pow(10, 3.5), pow(10, 4)]) * pow(10, 9)


    """
    get theta array and phi array ============================
    """

    phi = np.linspace(0, 2*pi, 70)

    print(np.degrees(phi))

    costheta = np.linspace(1, cos(radians(87)), 25)

    theta = np.arccos(costheta)

    print(np.degrees(theta))

    muon_v_dir = [[] for i in range(len(phi))]

    for i in range(0, len(theta)):
        for j in range(0, len(phi)):
            muon_v_dir[i].append(-1 * changecord(theta[i], phi[j]))



    for j in range(0, len(E)):
        print("energy: " + str(E[j]))
        p = open("2.4www5km/muon_2D_def_E_" + str(j) + "_pos.txt", "w")
        for i in range(0, len(muon_v_dir)):
            for k in range(0, len(muon_v_dir[i])):

                totalang1, totalang2, totaldis, energy, pos, prepos, predis, preenergy = main(muon_v_dir[i][k], p, E[j])

                p.write(str(theta[i]) + " " + str(phi[k]) + " " + str(totalang2) + " " + str(totaldis) + " " + str(energy) + " " + str(pos[0]) + " " + str(pos[1]) + " " + str(pos[2]) + " " + str(prepos[0]) + " " + str(prepos[1]) + " " + str(prepos[2]) + " " + str(predis) + " " + str(preenergy) + '\n')

                print("finish theta: " + str(degrees(theta[i])) + " phi: " + str(degrees(phi[k])) + "\n")
                print(str(theta[i]) + " " + str(phi[k]) + " " + str(totalang2) + " " + str(totaldis) + " " + str(energy) + " " + str(pos[0]) + " " + str(pos[1]) + " " + str(pos[2]) + " " + str(prepos[0]) + " " + str(prepos[1]) + " " + str(prepos[2]) + " " + str(predis) + " " + str(preenergy) + '\n')

        p.close()


start()

#depth()