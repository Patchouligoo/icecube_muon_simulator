import matplotlib.pyplot as plt
from math import *
from mpl_toolkits.mplot3d import Axes3D
import pylab
import numpy as np
from matplotlib import cm
import matplotlib.colors as colors
from scipy.interpolate import UnivariateSpline
np.set_printoptions(threshold=np.nan)



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


def getEloss(E):
    if E > pow(10, 28):
        print(E)
        print("!!!!!!!!!!!!!!!!! out of range")
        exit()

    E = E / pow(10, 9)

    # print(E)

    # print(spl_ionization(E))

    totalloss = (pow(10, spl_brems(log10(E))) * 2 + pow(10, spl_ionization(log10(E))) + pow(10, spl_photo(log10(E)))) * 0.9167 * pow(10, 9)

    return totalloss







def getratio(E):
    E = E/pow(10, 9)
    """
    ??????????????????????????????????????!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    """
    if E > 12 * pow(10, 4):
        return 1.4

    x = [97.82625165318657, 133.7675299099177, 200.77823789573986, 275.72398329251484, 352.1110743386039,
         476.5383449564929, 635.6351053049032,
         847.8271150190255, 1130.7166822094173, 1508.3206821038357, 2011.9989657768444, 2663.46265783033,
         4027.8127093079784, 5166.515006984569,
         7653.189919920697, 10210.089678867402, 13622.670782258743, 18175.964940003934, 24249.15815739906,
         32355.013124451427, 43166.44787435703,
         57592.638961903715, 76834.79570848323, 102514.57239865992, 119973.19012320979]
    y = [1.2611490711739954, 1.2677263867121407, 1.2750750472164714, 1.2804580767900025, 1.288172621958387,
         1.2983324829410647, 1.3073741724173016, 1.3172274171282663,
         1.3309291522102389, 1.3380521941185124, 1.3456334407705683, 1.3577839825843034, 1.3653796334056616,
         1.3677881117146828, 1.3778173204886548, 1.3816894036440746,
         1.3821837526278595, 1.3824920444548021, 1.3854757462098406, 1.385063818608925, 1.3876403768935077,
         1.3890337107410278, 1.3926461654031443, 1.3935212709862572, 1.3945520928716182]

    spl = UnivariateSpline(x, y, k= 4 ,s= 0.0001)

    return spl(E)

def getnewE_p(E, dis):

    deltax = dis / 20
    num = 20
    if dis > 1000:
        deltax = dis / 200
        num = 200

    if dis > 5000:
        deltax = dis / 500
        num = 500

    for i in range(0, num):
        E_loss1 = getEloss(E)
        E_temp = E + 100 * deltax * E_loss1
        E_loss2 = getEloss(E_temp)
        E_loss = (E_loss1 + E_loss2)/2

        E = E + 100 * deltax * E_loss
    return E


def getnewE_m(E, dis):
    deltax = dis / 20
    num = 20
    if dis > 1000:
        deltax = dis / 200
        num = 200

    if dis > 5000:
        deltax = dis / 500
        num = 500

    for i in range(0, num):
        E_loss1 = getEloss(E)
        E_temp = E - 100 * deltax * E_loss1
        E_loss2 = getEloss(E_temp)
        E_loss = (E_loss1 + E_loss2)/2

        E = E - 100 * deltax * E_loss


    return E
"""===================================================================================================================================================="""


lowerlimit = 0 #0
upperlimit = 25 #28

numtheta = 25 #25
numphi = 70


numberofthetainterest = numtheta

array_phi_theta_flux_total = np.array([])

arrayphi = []
arraytheta = []

for dis_num in range(0, 3):

    depth = -1

    if dis_num is 0:
        depth = 1.45
    if dis_num is 1:
        depth = 1.95
    if dis_num is 2:
        depth = 2.45

    for num in range(0, 6):
        check = 0

        if num < 3:
            check = 1
            f = open(str(depth) + "km/muon_2D_def_E_" + str(num) + "_pos.txt", "r")
        else:
            check = 2
            f = open(str(depth) + "km/muon_2D_def_E_" + str(num - 3) + "_neg.txt", "r")


        """create 2D list for def_ang, E and dis"""

        array_phi_theta_flux = [[] for k in range(numtheta)]

        """list for phi and theta"""

        for i in range(0, numtheta):
            arrayflux = []
            c = 0

            if i > numberofthetainterest:
                break
            else:
                for j in range(0, numphi):

                    temp = f.readline().split()

                    if i == 0 and num == 0 and dis_num == 0:
                        arrayphi.append(float(temp[1]))
                    if c == 0 and num == 0 and dis_num == 0:
                        arraytheta.append(cos(float(temp[0])))
                        c = c + 1

                    #E = float(temp[4])
                    dis = float(temp[3])

                    predis = float(temp[11])

                    pos = np.array([float(temp[5]), float(temp[6]), float(temp[7])])

                    prepos = np.array([float(temp[8]), float(temp[9]), float(temp[10])])

                    preenergy = float(temp[12])

                    """-----------------------------------------------"""
                    norm = np.linalg.norm(pos - prepos)
                    dir = (pos - prepos) / norm

                    #check_1 = 0
                    try:
                        step = norm / 2
                        pre_prepos = prepos
                        #print(E)

                        while abs(np.linalg.norm(prepos) - 6371000) > 0.01:
                            #print("in this iteration: -------------------------------------------")
                            #print("norm: " + str(norm) + " step: " + str(step))
                            #print(getEloss(preenergy))
                            #print(preenergy)

                            if np.linalg.norm(prepos) - 6371000 > 0:
                                #print(1)
                                prepos = prepos - dir * step
                                predis = predis - step
                                #e_temp = preenergy - step * getEloss(preenergy) * 100
                                #print(getEloss(E))
                                preenergy = getnewE_m(preenergy, step)
                                step = step / 2

                                #print(E_temp)
                                #print(getEloss(E_temp))
                                #print((getEloss(E) + getEloss(E_temp))/2 *step * 100)
                            else:
                                #print(2)
                                prepos = prepos + dir * step
                                predis = predis + step
                                #e_temp = preenergy + step * getEloss(preenergy) * 100
                                preenergy = getnewE_p(preenergy, step)
                                step = step / 2
                                #print(getEloss(E))
                                #print(E_temp)
                                #print(getEloss(E_temp))
                                #print((getEloss(E) + getEloss(E_temp))/2 *step * 100)
                        #e_temp = preenergy + getEloss(preenergy) * 100 * np.linalg.norm(prepos - pre_prepos)
                        #preenergy = preenergy + (getEloss(preenergy) + getEloss(e_temp)) / 2 * 100 * np.linalg.norm(prepos - pre_prepos)
                        print("theta: " + str(degrees(float(temp[0]))) + " E: " + str(preenergy) + " dis: " + str(predis))
                        # print(str(j) + " " + str(np.linalg.norm(pos) - 6371000))
                    except:
                        print(j)
                    """-----------------------------------------------"""
                    pos_normal = pos / np.linalg.norm(pos)

                    coss = np.dot(dir, pos_normal) / (np.linalg.norm(dir) * np.linalg.norm(pos_normal))

                    if check == 1:
                        #print("E is " + str(E/pow(10, 9)) + "GeV, and ratio is: " + str(getratio(E)))
                        Eflux = getratio(preenergy) * pow(preenergy / pow(10, 13), -3.7)
                    else:
                        Eflux = pow(preenergy / pow(10, 13), -3.7)
                    #if E > 1e+30: check_1 = 1

                    arrayflux.append(pow(coss, 2) * Eflux)

                """append np array to 2D list"""
                #avg = np.average(np.array(arrayflux))
                #if check_1 != 1:
                array_phi_theta_flux[i] = np.array(arrayflux)
                #else:
                    #exit()
                    #array_phi_theta_flux[i] = np.array([1 for i in range(numphi)])
        #print(np.array(array_phi_theta_flux))
        if num == 0 and dis_num == 0:
            array_phi_theta_flux_total = np.array(array_phi_theta_flux)
        else:
            array_phi_theta_flux_total = array_phi_theta_flux_total + np.array(array_phi_theta_flux)




arrayphi = np.array(arrayphi)
arraytheta = np.array(arraytheta)

#print(array_phi_theta_flux_total)

"""change 2D list to 2D np array"""
for i in range(0, len(array_phi_theta_flux_total)):
    #print(array_phi_theta_flux_total[i])
    #print(np.average(array_phi_theta_flux_total[i]))
    array_phi_theta_flux_total[i] = array_phi_theta_flux_total[i] / np.average(array_phi_theta_flux_total[i])

array_phi_theta_flux = array_phi_theta_flux_total


w = open("flux_data.txt", "w")

w.write(str(numtheta) + '\n')
w.write(str(numphi) + '\n')

for i in range(0, len(array_phi_theta_flux)):
    for j in range(0, len(array_phi_theta_flux[i])):
        w.write(str(array_phi_theta_flux[i][j]) + " ")

    w.write('\n')

for i in range(0, numtheta):
    w.write(str(arraytheta[i]) + " ")

w.write('\n')

for i in range(0, numphi):
    w.write(str(arrayphi[i]) + " ")

w.close()


