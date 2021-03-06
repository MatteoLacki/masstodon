import json
from collections import defaultdict


isotopes = [['Pr', [[140.90766200000002, 1.0]]],
            ['Ni',
             [[57.935342299999995, 0.680769095231328],
              [59.9307863, 0.262230419610671],
                 [60.931056299999995, 0.0113990830357779],
                 [61.928345400000005, 0.036346250253449],
                 [63.9279674, 0.0092551518687743]]],
            ['Yb',
             [[167.933892, 0.0012329299695777302],
              [169.934772, 0.0298222060986936],
                 [170.936332, 0.140905996539397],
                 [171.93639199999998, 0.21680068572105102],
                 [172.938222, 0.16102725365199302],
                 [173.938872, 0.32024990980512297],
                 [175.94258200000002, 0.129961018214165]]],
            ['Pd',
             [[101.905602, 0.0102075501879549],
              [103.9040311, 0.11146324882028301],
                 [104.90508090000002, 0.223336399264177],
                 [105.90348090000002, 0.27326441654003],
                 [107.9038929, 0.264546508837879],
                 [109.9051726, 0.11718187634967599]]],
            ['Pt',
             [[189.959934, 0.00012198734991181398],
              [191.961042, 0.007821588901230941],
                 [193.9626817, 0.32860592356572604],
                 [194.9647927, 0.33778897128367796],
                 [195.96495269999997, 0.25210785641529],
                 [197.967892, 0.0735536724841634]]],
            ['Ru',
             [[95.9075903, 0.0554029748080132],
              [97.90529599999999, 0.0187262734715792],
                 [98.9059348, 0.127588609866637],
                 [99.90421479999999, 0.126054915071901],
                 [100.90557790000001, 0.17058605337537802],
                 [101.9043449, 0.31545122520618396],
                 [103.90543199999999, 0.18618994820030801]]],
            ['Na', [[22.989769281999997, 1.0]]],
            ['Nb', [[92.90637199999999, 1.0]]],
            ['Nd',
             [[141.907732, 0.271519166958828],
              [142.909822, 0.121740433020292],
                 [143.910092, 0.23797766399758105],
                 [144.91258200000001, 0.0829297238509155],
                 [145.913122, 0.17189014035550199],
                 [147.916902, 0.0575610754128577],
                 [149.920902, 0.056381796404024]]],
            ['Mg',
             [[23.985041701, 0.789876809855212],
              [24.98583703, 0.100001999840013],
                 [25.98259302, 0.110121190304776]]],
            ['Li',
             [[6.0151228871, 0.0759339252859771], [7.016003443, 0.9240660747140229]]],
            ['Dy',
             [[155.924282, 0.000562985756460362],
              [157.924422, 0.00095297588970999],
                 [159.925202, 0.0232912107323685],
                 [160.926942, 0.18888942109764603],
                 [161.926812, 0.254747154896981],
                 [162.928742, 0.248957901365095],
                 [163.929182, 0.282598350261738]]],
            ['Y', [[88.905842, 1.0]]],
            ['Tl',
             [[202.97234509999998, 0.29520409591808194],
              [204.97442809999998, 0.7047959040819181]]],
            ['Tm', [[168.934222, 1.0]]],
            ['Rb',
             [[84.911789743, 0.7216911323547059],
              [86.90918053600001, 0.278308867645294]]],
            ['Ti',
             [[45.9526283, 0.0825200975882894],
              [46.9517593, 0.0744110706715194],
                 [47.947942299999994, 0.7371415430148379],
                 [48.947866299999994, 0.0541135063792345],
                 [49.944787299999994, 0.0518137823461185]]],
            ['Te',
             [[119.90406200000001, 0.0009097643710279039],
              [121.90304099999999, 0.0255053941029273],
                 [122.904271, 0.00892768772887822],
                 [123.902821, 0.047401722953755],
                 [124.904431, 0.0706966895574046],
                 [125.90331100000002, 0.18837621056146497],
                 [127.90446170000001, 0.317407791382032],
                 [129.906222759, 0.34077473934251]]],
            ['Rh', [[102.905502, 1.0]]],
            ['Ta',
             [[179.947462, 0.000120131992311553], [180.948002, 0.9998798680076879]]],
            ['Be', [[9.01218316, 1.0]]],
            ['Xe',
             [[123.905892, 0.000952296533640618],
              [125.904303, 0.0008901967596837949],
                 [127.9035318, 0.0191028304656971],
                 [128.904780864, 0.264005869018637],
                 [129.90350941, 0.0407099818156662],
                 [130.9050842, 0.21232352714236102],
                 [131.904155094, 0.269085350529324],
                 [133.90539569999999, 0.104356830141138],
                 [135.907214488, 0.08857311759385199]]],
            ['Pa', [[231.035882, 1.0]]],
            ['Ba',
             [[129.906322, 0.00106098514620795],
              [131.9050618, 0.00101098584619815],
                 [133.90450819999998, 0.0241714615995376],
                 [134.9056882, 0.0659202771161204],
                 [135.90457619999998, 0.0785413004217941],
                 [136.90582719999998, 0.112320827508415],
                 [137.9052472, 0.7169741623617271]]],
            ['Tb', [[158.925352, 1.0]]],
            ['La',
             [[137.90712299999998, 0.0008881718721032509],
              [138.906362, 0.9991118281278969]]],
            ['Si',
             [[27.976926535300002, 0.9222208333499999],
              [28.976494665300002, 0.0468584376987476],
                 [29.973770012, 0.0309207289512526]]],
            ['As', [[74.92159570000001, 1.0]]],
            ['U',
             [[234.040952, 5.4599923560107e-05],
              [235.04393199999998, 0.00720468991343412],
                 [238.050792, 0.992740710163006]]],
            ['W',
             [[179.94671200000002, 0.00120987296333885],
              [181.9482047, 0.264988176241495],
                 [182.9502237, 0.143124971877953],
                 [183.95093169999998, 0.306387829277926],
                 [185.954362, 0.28428914963928803]]],
            ['Gd',
             [[151.919802, 0.00200963625583769],
              [153.920872, 0.0218260494850432],
                 [154.922632, 0.147985214676144],
                 [155.922132, 0.20467295419529102],
                 [156.923972, 0.156491675006824],
                 [157.924112, 0.24843503325898003],
                 [159.927062, 0.21857943712188102]]],
            ['Fe',
             [[53.9396093, 0.0584527927212081],
              [55.934936300000004, 0.9175324978567758],
                 [56.935393299999994, 0.0211907435920025],
                 [57.9332743, 0.00282396583001346]]],
            ['Br', [[78.9183381, 0.5068988961766121],
                    [80.9162901, 0.4931011038233879]]],
            ['Sr',
             [[83.91341990000001, 0.00560977560897564],
              [85.9092619, 0.0986060557577697],
                 [86.9088789, 0.0700071997120115],
                 [87.9056139, 0.825776968921243]]],
            ['Hf',
             [[173.940052, 0.0016096523150999401],
              [175.941412, 0.0526686235773073],
                 [176.943232, 0.18596983051660798],
                 [177.943712, 0.27282107064874],
                 [178.945822, 0.136190582834108],
                 [179.946562, 0.350740240108137]]],
            ['Mo',
             [[91.90680859999999, 0.145308494342837],
              [93.9050853, 0.0914964585241384],
                 [94.9058393, 0.15838755864132098],
                 [95.90467629999999, 0.166690329831185],
                 [96.9060183, 0.0959997920307794],
                 [97.9054053, 0.243900902666405],
                 [99.9074728, 0.0982164639633334]]],
            ['He',
             [[3.016029322, 1.342999991942e-06], [4.00260325414, 0.9999986570000079]]],
            ['C', [[12.0, 0.989211941850467], [13.0033548352, 0.0107880581495331]]],
            ['B',
             [[10.0129373, 0.199480830670927], [11.009305300000001, 0.8005191693290741]]],
            ['F', [[18.998403163699997, 1.0]]],
            ['I', [[126.904473, 1.0]]],
            ['H',
             [[1.0078250322700002, 0.9998842901643079],
              [2.01410177819, 0.000115709835692033]]],
            ['K',
             [[38.963706493000004, 0.9325805260710841],
              [39.96399824, 0.00011709988524211201],
                 [40.961825263, 0.0673023740436734]]],
            ['Mn', [[54.938044299999994, 1.0]]],
            ['Lu', [[174.940782, 0.9740087675772041],
                    [175.942692, 0.0259912324227957]]],
            ['O',
             [[15.9949146202, 0.9975676097295609],
              [16.9991317576, 0.000380998476006096],
                 [17.999159613699998, 0.00205139179443282]]],
            ['Ne',
             [[19.992440182, 0.9047666663333571],
              [20.99384673, 0.0027098103132780697],
                 [21.99138512, 0.0925235233533653]]],
            ['P', [[30.9737619986, 1.0]]],
            ['S',
             [[31.972071174099998, 0.9498500119990401],
              [32.971458910100004, 0.00751939844812415],
                 [33.96786703, 0.0425205983521318],
                 [35.9670812, 0.000109991200703944]]],
            ['Th', [[232.03806200000002, 1.0]]],
            ['Re',
             [[184.9529559, 0.37400503979840793], [186.955751, 0.6259949602015921]]],
            ['Kr',
             [[77.9203656, 0.00355294812695735],
              [79.9163786, 0.022860666234273],
                 [81.9134837, 0.115931407401452],
                 [82.9141272, 0.11500022099677303],
                 [83.911497733, 0.569863179997572],
                 [85.91061063299999, 0.172791577242972]]],
            ['Sm',
             [[143.912012, 0.0307725222770867],
              [146.914902, 0.149881578776357],
                 [147.91483200000002, 0.112382691006086],
                 [148.917192, 0.138246406123312],
                 [149.917282, 0.0737920685273479],
                 [151.91974199999999, 0.267451009404715],
                 [153.922222, 0.227473723885096]]],
            ['V',
             [[49.9471567, 0.0025039799681602602],
              [50.943957700000006, 0.9974960200318399]]],
            ['Sc', [[44.9559086, 1.0]]],
            ['Sb',
             [[120.90381200000002, 0.572091349038115],
              [122.90421200000002, 0.427908650961885]]],
            ['Bi', [[208.980401, 1.0]]],
            ['N',
             [[14.0030740042, 0.996358014567942], [15.0001088994, 0.00364198543205827]]],
            ['Os',
             [[183.95248909999998, 0.00020994772301696898],
              [185.95384099999998, 0.0159260344174301],
                 [186.955751, 0.0196151158361568],
                 [187.955841, 0.132457018202468],
                 [188.95814199999998, 0.161519781574388],
                 [189.958442, 0.262554623898649],
                 [191.96148200000002, 0.40771747834789096]]],
            ['Se',
             [[73.92247591, 0.00893842683687671],
              [75.91921372, 0.0937125065988386],
                 [76.91991426, 0.0763025707475484],
                 [77.9173092, 0.23768616723456698],
                 [79.9165229, 0.49605369454975895],
                 [81.9167001, 0.0873066340324103]]],
            ['Hg',
             [[195.965832, 0.0015098158024721],
              [197.9667693, 0.0997078356440514],
                 [198.9682813, 0.168701418426952],
                 [199.9683273, 0.230990819120067],
                 [200.9703036, 0.131793921141621],
                 [201.9706436, 0.29858957207220693],
                 [203.97349430000003, 0.0687066177926293]]],
            ['Zn',
             [[63.9291426, 0.49164571388582],
              [65.9260347, 0.27732550874018397],
                 [66.9271287, 0.0404052925974617],
                 [67.9248457, 0.184515103497573],
                 [69.925322, 0.0061083812789610795]]],
            ['Co', [[58.933194400000005, 1.0]]],
            ['Ag', [[106.905092, 0.518389668985958],
                    [108.9047551, 0.48161033101404205]]],
            ['Cl', [[34.96885273, 0.7575948481030379],
                    [36.96590264, 0.242405151896962]]],
            ['Ca',
             [[39.96259092, 0.9694008384267271],
              [41.9586181, 0.00647222841715371],
                 [42.9587662, 0.00135098505810526],
                 [43.955482200000006, 0.0208608692787858],
                 [45.953692, 4.2999524425259796e-05],
                 [47.95252289, 0.001872079294803]]],
            ['Ir', [[190.960592, 0.37305077968812506],
                    [192.962922, 0.6269492203118749]]],
            ['Eu',
             [[150.91986200000002, 0.47810306557082], [152.921242, 0.5218969344291801]]],
            ['Al', [[26.98153858, 1.0]]],
            ['Ce',
             [[135.9071293, 0.00185197333158403],
              [137.905998, 0.00251196382772088],
                 [139.905442, 0.8844924633085279],
                 [141.909252, 0.111143599532167]]],
            ['Cd',
             [[105.90646090000001, 0.0125671975149542],
              [107.90418390000002, 0.00892800905398096],
                 [109.9030074, 0.12489014949666198],
                 [110.90418340000001, 0.12798345968848898],
                 [111.90276340000001, 0.24126719741497601],
                 [112.90440829999999, 0.12218475280012599],
                 [113.90336529999999, 0.287277937020045],
                 [115.9047632, 0.0749012970107666]]],
            ['Ho', [[164.930332, 1.0]]],
            ['Ge',
             [[69.9242497, 0.20570581230133303],
              [71.92207586, 0.27450372611621],
                 [72.92345904, 0.0775040170862401],
                 [73.92117776100001, 0.3649824068120979],
                 [75.92140272, 0.0773040376841185]]],
            ['Ar',
             [[35.96754512, 0.0033362057963807],
              [37.962732200000005, 0.0006297992064530001],
                 [39.962383122, 0.9960339949971659]]],
            ['Au', [[196.96656959999999, 1.0]]],
            ['Zr',
             [[89.904702, 0.51442271162175],
              [90.905642, 0.112234410554394],
                 [91.90503199999999, 0.171550886397901],
                 [93.906312, 0.17378837625021498],
                 [95.908272, 0.0280036151757399]]],
            ['Ga', [[68.9255749, 0.6010797978404039], [70.9247037, 0.398920202159596]]],
            ['In',
             [[112.90406270000001, 0.0429548454185498],
              [114.90387878899999, 0.9570451545814499]]],
            ['Cs', [[132.905451967, 1.0]]],
            ['Cr',
             [[49.9460427, 0.04345074383047899],
              [51.940506400000004, 0.837881075122238],
                 [52.94064839999999, 0.0950104838658065],
                 [53.9388794, 0.0236576971814761]]],
            ['Pb',
             [[203.9730449, 0.014094362255098001],
              [205.9744669, 0.24100359856057602],
                 [206.9758979, 0.221011595361855],
                 [207.9766539, 0.523890443822471]]],
            ['Er',
             [[161.928792, 0.00139597347650395],
              [163.929212, 0.0160126957587806],
                 [165.930302, 0.335027234482545],
                 [166.932052, 0.22868665495355603],
                 [167.93238200000002, 0.26977667424318896],
                 [169.935472, 0.149100767085425]]],
            ['Cu',
             [[62.929598399999996, 0.6914942551723451], [64.9277906, 0.308505744827655]]],
            ['Sn',
             [[111.9048244, 0.00970737900766793],
              [113.9027837, 0.006608215781738929],
                 [114.90334471, 0.0034090795485219],
                 [115.90174309999999, 0.14537074989752802],
                 [116.9029543, 0.0768592480030392],
                 [117.90160729999998, 0.24214462095234301],
                 [118.9033116, 0.0859168024633349],
                 [119.9022027, 0.325722055045138],
                 [121.903442, 0.0463174942765453],
                 [123.9052778, 0.0579443550241435]]]]


def get_isotopic_masses_and_probabilities():
    """Retrieve the information on masses and frequencies of isotopes.

    Returns
    =======
    out : tuple
        Two dictionaries, both with keys set to element encodings,
        such as Ag, Al, Ar, ...
        The values of the first dictionary are lists of masses of elements.
        The values of the first dictionary are lists of probabilities of elements.
    """
    iso_masses = defaultdict(list)
    iso_probs  = defaultdict(list)
    for element, isos in isotopes:
        for mass, prob in isos:
            iso_masses[element].append(mass)
            iso_probs[element].append(prob)
    return iso_masses, iso_probs


def get_isotopic_masses():
    """Retrieve the information on masses of isotopes."""
    iso_masses = defaultdict(list)
    for element, isos in isotopes:
        for mass, prob in isos:
            iso_masses[element].append(mass)
    return iso_masses


def get_isotopic_probabilities():
    """Retrieve the information on frequencies of isotopes."""
    iso_probs  = defaultdict(list)
    for element, isos in isotopes:
        for mass, prob in isos:
            iso_probs[element].append(prob)
    return iso_probs

