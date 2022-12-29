logpi = 44
w = 2^logpi
sigma=3.2
alpha=36


# lwe - RLWE, MLWE
# t - plaintext size
# n - ring size for MLWE or hint on ring size for RLWE
# d - hint on MLWE rank not used on RLWE
confs = (
         {'s':'RLWE', 't':2   ,  'n':16384, 'd':1},
         {'s':'MLWE', 't':2   ,   'n':128, 'd':89},
         {'s':'RLWE', 't':2^16, 'n':32768, 'd':1},
         {'s':'MLWE', 't':2^16,   'n':128, 'd':180},
         {'s':'MLWE', 't':2^16,   'n':512, 'd':45},
        )


Berr = sigma * sqrt(alpha)


def delta(n):
    return (2. * sqrt(n))
def Vnorm(n):
    return Berr * (1. + 2. * delta(n) * Bkey)
def epsilon1(n):
    return 5 / (delta(n) * Bkey)
def C1(n):   
    return (1 + epsilon1(n)) * delta(n) * delta(n) * t * Bkey
def C2(n,loq):
    return delta(n) * delta(n) * Bkey * ((1 + 0.5) * Bkey + t * t) + delta(n) * (floor(logq / (log(2) * logpi)) + 1) * w * Berr
def findLogq(n, logq):
    #print("n: ",n,"C1:",RR(C1(n)),"C2:",RR(C2(n,logq)), "t: ",t,"l:",l,"Vnorm:",RR(Vnorm(n)))
    #print("delta: ", RR(delta(n)), "epsilon: ",RR(epsilon1(n)))
    return log(4 * t) + (l - 1) * log(C1(n)) + log(C1(n) * Vnorm(n) + l * C2(n, logq))

def findMLWEdelta(n, logq):
    #print("findMLWEdelta",n,RR(logq))
    load("https://bitbucket.org/malb/lwe-estimator/raw/HEAD/estimator.py")
    q = 2^ceil(logq/log(2))
    #set_verbose(0)
    Logging.set_level(Logging.NOTSET)
    alpha = alphaf(sigmaf(sigma),q)
    L = estimate_lwe(n, alpha, q, reduction_cost_model=BKZ.enum)
    delta_enum1 = L['usvp']['delta_0'] 
    delta_enum2 = L['dec']['delta_0']  
    delta_enum3 = L['dual']['delta_0']  
    L = estimate_lwe(n, alpha, q, reduction_cost_model=BKZ.sieve)
    delta_sieve1 = L['usvp']['delta_0'] 
    delta_sieve2 = L['dec']['delta_0']  
    delta_sieve3 = L['dual']['delta_0']
    return max(delta_enum1,delta_enum2,delta_enum3,delta_sieve1,delta_sieve2,delta_sieve3)

def findMLWEn(mlwe, d, logq):               # dimension of the Module-LWE problem
    #print("findMLWEn",mlwe,d,RR(logq))
    mlwe_hardness = findMLWEdelta(mlwe*d,logq)
    while mlwe_hardness > 1.0045:           # increasing the mlwe dimension until MLWE provides ~ 128-bit security
        if(schema=='RLWE'):
            mlwe *= 2
        else:
            d += 1
        mlwe_hardness = findMLWEdelta(mlwe*d, logq)
    return mlwe, d

def findLogQ(conf, logq):
    n = conf['n']
    d = conf['d']
    #print("findLogQ",n,d,RR(logq))
    logq = findLogq(n*d, logq)
    nn,nd = findMLWEn(n, d, logq)
    while (nn > n or nd > d):
        n = nn
        d = nd
        logq = findLogq(n*d, logq)
        nn,nd = findMLWEn(n, d, logq)
    logq_new = findLogq(n*d, logq)
    while (abs(logq_new - logq) > log(1.001)):
        #print("Last",n,d,RR(logq))
        logq = logq_new
        logq = findLogq(n*d, logq)
    return conf['t'],n,d,logq/log(2)+1

results = []
for conf in confs:
    l = max(log(conf['t'])/log(2)+1,8)
    logq = 6. * log(10);
    schema = conf['s']
    if(schema=='RLWE'):
        Bkey = sigma * sqrt(alpha)
    else:
        Bkey=1
    results.append(findLogQ(conf, logq))

ref_space = {2:128,65536:256}
ref_speed = {2:16384,65536:32768}
print("    t|","    n|","  d|","logq|" ,"ncrt|","  k|"," s. imp|","slowdown")
print("-----+------+----+-----+-----+----+--------+---------")
for r in results:
    space_imp = (ref_space[r[0]]-(r[1]*r[2])/128)/ref_space[r[0]]*100
    speed_slowdown = (r[2]*r[2]*r[1])/ref_speed[r[0]]
    print("%5d| %5d| %3d| %4d| %4d| %3d| %7.2f| %5d" % (r[0],r[1],r[2],ceil(r[3]),ceil(r[3]/logpi),ceil(r[1]/128),space_imp,ceil(speed_slowdown)))

