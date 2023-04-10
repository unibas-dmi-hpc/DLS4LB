from mpi4py import MPI
import random
import math

WRK_TAG = 9970
REQ_TAG = 9960
END_TAG = 9940


DLS_Schedulers = {
    0:"STATIC",
    1:"SS",
    2:"FSC",
    3:"mFSC",
    4:"GSS",
    5:"TSS",
    6:"FAC",
    7:"WF",
    8:"AWF",
    9:"AWF_B",
    10:"AWF_C",
    11:"AWF_D",
    12:"AWF_E",
    13:"AF",
    14:"TFSS",
    15:"FISS",
    16:"VISS",
    17:"RND",
    18:"PLS",
    19:"TAP"
}

class infoDLS:
    nProcs = None
    myRank = None
    foreman = None
    firstIter = None
    lastIter = None
    method = None
    chunkStart = None
    wStart = None
    itersScheduled = None          
    batchSize = None
    batchRem = None
    numENDed = None
    minChunkSize = None
    maxChunkSize = None
    numENDed = 0
    t0 = t1 = 0
    sumt1 = sumt2 = 0
    rStart = rSize = 0
    def __init__(self, comm, size, rank, firstIter, lastIter, master, method, requestWhen, breakAfter, minChunk, h_overhead, sigma, mu, alpha, nKNL, Xeon_speed, KNL_speed, X, B, SWR):
        # self.comm = MPI.COMM_WORLD
        self.comm = comm
        self.nProcs = size
        self.myRank = rank
        # self.nProcs = self.comm.Get_size()
        # self.myRank = self.comm.Get_rank()
        self.firstIter = firstIter
        self.lastIter = lastIter
        self.N = self.lastIter - self.firstIter + 1
        self.foreman = master
        self.firstRank = 0
        self.lastRank = size - 1
        # self.lastRank = self.nProcs - 1
        self.wSize = 0
        self.gotWork = 1

        self.method = method
        if method > len(DLS_Schedulers) and method < 0:
            self.method = 0

        self.stats = [0] * (3 * size)
        self.weights = [0] * size
        self.requestWhen = requestWhen
        self.breakAfter = breakAfter
        self.workTime = 0.0
        self.myIters = 0
        self.timeStep = 1
        self.nextWRKrcvd = 0
        self.req4WRKsent = 0
        self.finishedOne = 0
        self.probeFreq = max(1, breakAfter)
        self.sendRequest = max(1, requestWhen)
        self.minChunk = minChunk
        self.h_overhead = h_overhead
        self.sigma = sigma
        self.numChunks = 0
        self.subChunkSize = 0
        self.myExecs = 0
        self.mySumTimes = 0.0
        self.mySumSizes = 0.0
        totalSum = nKNL * KNL_speed + (self.nProcs - nKNL) * Xeon_speed

        for i in range(0, size):
            self.stats[3*i] = -1    # mu
            self.stats[3*i+1] = -1  # sigma
            self.stats[3*i+2] = 0   # performance data count

            # initializing weights
            if i < (size - nKNL):
                coreSpeed = Xeon_speed
            else:
                coreSpeed = KNL_speed

            self.weights[i] = coreSpeed/totalSum * self.nProcs

        # FSC inital chunk calc.
        k_FSC_num = (math.sqrt(2) * self.N * self.h_overhead) 
        k_FSC_denom =  (self.sigma * self.nProcs * math.sqrt(math.log(self.nProcs)))
        self.chunkFSC = 0
        if k_FSC_denom != 0:
            self.chunkFSC = k_FSC_num / k_FSC_denom
        if self.chunkFSC == 0:
            self.chunkFSC = minChunk
        self.chunkFSC = math.pow(self.chunkFSC, (2.0/3.0))
        self.chunkFSC = int(math.ceil(self.chunkFSC))
        
        # MFSC inital chunk calc.
        MFSC_size = (self.N + self.nProcs - 1) / self.nProcs
        self.chunkMFSC = math.floor((0.55 + MFSC_size * math.log(2)) / math.log(MFSC_size))

        # AWF
        if self.method == 8 and self.myRank == self.foreman:
            if self.timeStep == 1:
                for i in range(self.firstRank, self.lastRank+1):
                    self.weights[i] = 1.0
            else:
                awap = 0.0
                for i in range(self.firstRank, self.lastRank+1):
                    awap += self.stats[3*i]
                
                awap = awap/self.nProcs
                trw = 0

                for i in range(self.firstRank, self.lastRank+1):
                    trw += awap/self.stats[3*i]

                for i in range(self.firstRank, self.lastRank+1):
                    self.weights[i] = ((awap/self.stats[3*i]) * self.nProcs) / trw

        # TSS
        self.chunkTSS = int(math.ceil(self.N / (2 * self.nProcs)))
        S = math.ceil((2 * self.N) / (self.chunkTSS + 1))
        self.deltaTSS = int(math.floor((self.chunkTSS - 1) / (S - 1)))

        # FISS
        self.chunkFISS = self.N / ((2 + B) * self.nProcs)
        fiss_num = 2 * self.N * (1 - (B / (2 + B)))
        fiss_denom = self.nProcs * B * (B - 1)
        self.deltaFISS = int(math.ceil(fiss_num / fiss_denom))

        # VISS
        self.VISSremain = math.floor(self.chunkFISS)
        self.chunkVISS = 0

        # PLS
        self.chunkPLS = (self.N * SWR) / self.nProcs
        self.PLSremain = self.N - (self.N * SWR)

        # TAP
        self.V_alpha = (alpha * sigma) / mu
        

def StartLoop(info):
    NULLloc = 0
    if info.myRank == info.foreman:
        info.chunkStart = info.firstIter
        info.itersScheduled = 0          
        info.batchSize = 0
        info.batchRem = 0
        info.numENDed = 0
        info.numChunks = 0
        
        if info.minChunk > 0:
            info.minChunkSize = info.minChunk
        else:
            info.minChunkSize = max(1, info.chunkMFSC / 2)
        
        info.maxChunkSize = (info.N + 2*info.nProcs) / (2*info.nProcs)

        for i in range(0, info.nProcs):    # i is a worker
            if info.chunkStart < info.lastIter:
                SendChunk(info, i)
            else:
                # MPI.COMM_WORLD.Send()
                info.comm.send(NULLloc, i, tag=END_TAG)
                info.numENDed += 1



def GetChunkSize(info, rank):
    
    rem = info.N - info.itersScheduled

    if info.method == 0: # STATIC
        tCh = info.batchSize = math.ceil( info.N / info.nProcs)
        info.batchRem = min(info.batchSize, rem)

    elif info.method == 1: # ss
        tCh = info.batchSize = 1
        info.batchRem = min(info.batchSize, rem)
    
    elif info.method == 2: # FSC
        tCh = info.batchSize = min(info.chunkFSC, rem)
        info.batchRem = min(info.batchSize, rem)
    
    elif info.method == 3: # mFSC
        tCh = info.batchSize = min(info.chunkMFSC, rem)
        info.batchRem = min(info.batchSize, rem)
    
    elif info.method == 4: # GSS
        tCh = info.batchSize = min(rem, max((rem + info.nProcs - 1) // info.nProcs, info.minChunkSize))
        info.batchRem = min(info.batchSize, rem)

    elif info.method == 5: # TSS
        tCh = info.batchSize = max(info.minChunkSize, min(rem, info.chunkTSS))
        info.batchRem = min(info.batchSize, rem)
        info.chunkTSS = max(info.minChunkSize, min(rem, info.chunkTSS)) - info.deltaTSS

    elif info.method == 6: # FAC
        
        if info.batchRem == 0:
            info.batchSize = info.nProcs * max(info.minChunkSize, (rem + 2 * info.nProcs - 1)//(2 * info.nProcs))
            info.batchRem = min(info.batchSize, rem)
        
        tCh = min(rem, max(info.minChunkSize, info.batchSize // info.nProcs))
        
    elif info.method == 8 or info.method == 7: # AWF or WF
        
        if info.batchRem == 0:
            info.batchSize = info.nProcs * max(info.minChunkSize, (rem + 2 * info.nProcs - 1)//(2 * info.nProcs))
            info.batchRem = min(info.batchSize, rem)
        
        tCh = min(rem, max(info.minChunkSize, info.batchSize // info.nProcs*info.weights[rank]))
    
    elif info.method == 9 or info.method == 11: # AWF_B OR AWF_D
        tCh = 0
        if info.stats[3*rank] < 0:
            info.batchRem = info.batchSize = min(rem, info.minChunkSize)
        else: # all ranks have wap
            awap = 0.0
            for i in range(info.firstRank, info.lastRank+1):
                awap += info.stats[3*i]
            
            awap /= info.nProcs

            trw = 0.0 # total ref weight (refwt(i) = awap/info->stats[3*i] 
            for i in range(info.firstRank, info.lastRank+1):
                trw += (awap/info.stats[3*i]) 

            weight = ((awap/info.stats[3*rank]) * info.nProcs) / trw
            
            if info.batchRem == 0:
                tCh = max(info.minChunkSize, (rem+2*info.nProcs-1) / (2*info.nProcs))
                info.batchSize = info.nProcs * tCh
                info.batchRem = min(info.batchSize, rem)

            tCh = weight * (info.batchSize / info.nProcs) + 0.55
            tCh = max(info.minChunkSize, tCh)
            tCh = min(rem, tCh)
            
    elif info.method == 10 or info.method == 12: # AWF_C OR AWF_E
        tCh = 0
        if info.stats[3*rank] < 0.0:
            tCh = info.minChunkSize
        else: # all ranks have wap
            awap = 0.0
            for i in range(info.firstRank, info.lastRank+1):
                awap += info.stats[3*i]
            
            awap /= info.nProcs

            trw = 0.0 # total ref weight (refwt(i) = awap/info->stats[3*i] 
            for i in range(info.firstRank, info.lastRank+1):
                trw += (awap/info.stats[3*i]) 
            
            weight = ((awap/info.stats[3*rank]) * info.nProcs) / trw
            tCh = weight * ((rem+2*info.nProcs-1) / (2*info.nProcs)) + 0.55
            
        tCh = info.batchSize = max(info.minChunkSize, tCh)
        info.batchRem = min(rem, tCh)
        
    elif info.method == 13: # AF

        if info.stats[3*rank] < 0.0:
            tCh = info.minChunkSize

        else:
            bigD = 0.0
            bigT = 0.0

            for i in range(info.firstRank, info.lastRank+1):
                bigD += (info.stats[3*i+1]/info.stats[3*i])
                bigT += (1.0/info.stats[3*i])
            
            bigT = 1.0/bigT

            tCh = 0.55 + (0.5 * (bigD + 2.0 * bigT * rem - math.sqrt(bigD * (bigD + 4.0 * bigT * rem))) / info.stats[3*rank])
            tCh = min(info.maxChunkSize, tCh)

        tCh = info.batchSize = max(info.minChunkSize, tCh)
        info.batchRem = min(info.batchSize, rem)

    elif info.method == 14: # TFSS
        tCh = 0
        totalTSS = 0
        tempCh = 0
        # print(info.deltaTSS, " ", info.chunkTSS)
        if info.batchRem == 0:
            tCh = max(info.minChunkSize, info.chunkTSS)
            totalTSS += tCh
            tempCh = tCh
            # print(totalTSS, " ", tempCh)
            for i in range(1, info.nProcs):
                totalTSS += tempCh - info.deltaTSS
                tempCh -= info.deltaTSS
            info.chunkTSS = tCh - (info.nProcs * info.deltaTSS)
            if totalTSS > 0:
                info.batchSize = info.nProcs * int(math.floor(totalTSS / info.nProcs))
            else:
                info.batchSize = info.nProcs * info.minChunkSize
            
            info.batchRem = min(info.batchSize, rem)
        
        tCh = max(info.minChunkSize, int(info.batchSize / info.nProcs))
        tCh = min(rem, tCh)

    elif info.method == 15: # FISS
        tCh = 0
        if info.batchRem == 0:
            tCh = info.chunkFISS
            info.chunkFISS += info.deltaFISS
            tCh = min(rem, tCh)
            info.batchSize = tCh * info.nProcs
            info.batchRem = min(info.batchSize, rem)
        
        tCh = info.batchSize / info.nProcs
        tCh = int(min(rem, tCh))
    
    elif info.method == 16: # VISS
        tCh = 0
        if info.batchRem == 0:
            tCh = math.ceil(info.chunkVISS // 2) + info.VISSremain
            info.chunkVISS = tCh
            info.batchSize = tCh * info.nProcs
            info.batchRem = min(info.batchSize, rem)
        else:
            tCh = info.chunkVISS
        tCh = max(info.minChunkSize, tCh)
        tCh = min(rem, tCh)
    
    elif info.method == 17: # RND
        tCh = random.randint(1, 250)
        tCh = min(tCh, rem)
        info.batchSize = tCh
        info.batchRem = min(info.batchSize, rem)
    
    elif info.method == 18: # PLS
        if info.itersScheduled < info.N - info.PLSremain:
            tCh= info.chunkPLS
        else:
            tCh = max((rem + info.nProcs-1) / info.nProcs, info.minChunkSize )
            
        tCh = min (rem, tCh)
        info.batchSize = tCh
        info.batchRem = min(info.batchSize, rem)
    
    elif info.method == 19: # TAP
        tCh = max((rem + info.nProcs - 1) // info.nProcs, info.minChunkSize)
        V_alpha_2 = math.pow(info.V_alpha, 2)
        tCh += (V_alpha_2/2) - info.V_alpha * math.sqrt(2 * tCh + (V_alpha_2 / 4))
        info.batchSize = tCh
        info.batchRem = min(info.batchSize, rem)
    else:
        if info.method < 0 or info.method > 19:
            tCh = (info.N + info.nProcs - 1) / info.nProcs
            i = info.N % info.nProcs

            if i > 0 and rank >= 1:
                tCh -= 1
            
            tCh = min(tCh, rem)
            info.batchSize = tCh
            info.batchRem = min(info.batchSize, rem)
        
    # print(info.batchRem, " ", info.batchSize)
    chunksize = min(info.batchRem, tCh)

    info.batchRem = info.batchRem - chunksize

    if info.batchRem > 0 and info.batchRem <= info.minChunkSize:
        chunksize += info.batchRem
        info.batchRem = 0
    
    return int(chunksize)

def SendChunk(info, worker):

    chunksize = GetChunkSize(info, worker)
    chunkInfo = [0] * (2)
    # print(type(chunksize))
    # print(chunksize)
    chunkInfo[0] = info.chunkStart
    chunkInfo[1] = chunksize

    if worker == info.foreman:

        if info.wSize == 0:

            info.t0 = MPI.Wtime()
            info.wStart = chunkInfo[0]
            info.wSize = chunkInfo[1]
            info.rStart = info.wStart
            info.rSize = info.wSize
            info.req4WRKsent = 0

            SetBreaks(info)

            info.sumt1 = 0.0
            info.sumt2 = 0.0
        
        else:

            info.nextStart = chunkInfo[0]
            info.nextSize = chunkInfo[1]
            info.nextWRKrcvd = 1
    
    else:
        # print("sneding chunk information to ", worker, " from ", info.myRank)
        info.comm.send(chunkInfo, worker, tag=WRK_TAG)
    
    info.chunkStart += chunksize
    info.itersScheduled += chunksize

    info.numChunks += 1

def  SetBreaks (info):

    if (info.myRank == info.foreman):
        # /* when to check for messages */
        if info.breakAfter < 0:
            info.probeFreq = max( 1, (info.wSize+info.nProcs-1)//info.nProcs//4 )
        else:
            info.probeFreq = max( 1, info.breakAfter)

        # /* how many iterates left before requesting next chunk */
        if info.requestWhen < 0:
            info.sendRequest = info.probeFreq
        else:
            info.sendRequest = info.requestWhen
    
    else: #{ /* not the foreman */

        # /* how many iterates left before requesting next chunk */
        if info.requestWhen < 0:
            info.sendRequest = max( 1, (15*info.wSize)//100 )
        else:
            info.sendRequest = info.requestWhen

        # /* when to check for messages */
        info.probeFreq = max(1, info.wSize - info.sendRequest)


def DLS_Terminated(info):
    done = 0

    if info.N <= 0:
        done = 1
    else:
        done = (info.gotWork == 0) and (info.wSize == 0)
    return done

def DLS_StartChunk(info):

    chunkStart = chunkSize = 0
    NULLloc = None
    msgInQueue = 0
    loc = 0
    updateJ = 0
    chunkInfo = [0] * 2
    perfInfo = [0] * 4
    STAT = tStatus = MPI.Status()
    if info.comm != MPI.COMM_NULL:

        if info.wSize == 0:
            info.comm.probe(MPI.ANY_SOURCE, MPI.ANY_TAG, STAT)
            msgInQueue = 1
        else:
            msgInQueue = info.comm.iprobe(MPI.ANY_SOURCE, MPI.ANY_TAG, STAT)

        while msgInQueue:

            if STAT.Get_tag() == WRK_TAG:
  
                chunkInfo = info.comm.recv(source=STAT.Get_source(), tag=WRK_TAG, status=tStatus)
                if info.wSize == 0:
                    info.t0 = MPI.Wtime()
                    info.wStart = chunkInfo[0]
                    info.wSize = chunkInfo[1]
                    info.rStart = info.wStart
                    info.rSize = info.wSize
                    info.req4WRKsent = 0

                    SetBreaks(info)

                    info.sumt1 = 0.0
                    info.sumt2 = 0.0
        
                else:
                    info.nextStart = chunkInfo[0]
                    info.nextSize = chunkInfo[1]
                    info.nextWRKrcvd = 1

            elif STAT.Get_tag() == REQ_TAG:        
                worker = STAT.Get_source()
                perfInfo = info.comm.recv(source=worker, tag=REQ_TAG, status=tStatus)

                if info.method >= 9 and info.method <=13:
                    loc = int(perfInfo[2])
                    info.stats[3*loc+2] = info.stats[3*loc+2] + 1.0
                    info.stats[3*loc] = perfInfo[0]
                    info.stats[3*loc+1] = perfInfo[1]

                    if info.finishedOne != info.nProcs:
                        updateJ = loc
                        for i in range(info.firstRank, info.lastRank + 1):
                            if info.stats[3 * i + 2] > 0.0 and info.stats[3 * i] < info.stats[3 * updateJ]:
                                updateJ = i
                        
                        info.finishedOne = 0

                        for i in range(info.firstRank, info.lastRank + 1):

                            if info.stats[3*i+2] == 0.0:
                                info.stats[3*i] = info.stats[3*updateJ]
                                info.stats[3*i+1] = info.stats[3*updateJ+1]
                            
                            else:
                                info.finishedOne += 1

                if info.chunkStart <= info.lastIter:
                    SendChunk(info, worker)
                else:
                    info.numENDed += 1
                    if worker != info.myRank:
                        info.comm.send(NULLloc, worker, tag= END_TAG)

                    info.gotWork = info.numENDed != info.nProcs
            
            elif STAT.Get_tag() == END_TAG:
                NULLloc = info.comm.recv(source=STAT.Get_source(), tag=STAT.Get_tag(), status = tStatus)
                info.gotWork = 0
        
            msgInQueue = info.comm.iprobe(source=MPI.ANY_SOURCE, tag = MPI.ANY_TAG, status = STAT)
            
    
    chunkStart = info.wStart
    chunkSize = min(info.wSize, info.probeFreq)
    if info.method == 13:
        chunkSize = min(1, chunkSize)
    
    info.subChunkSize = chunkSize
    if info.subChunkSize != 0:
        info.t1 = MPI.Wtime()
    
    # print(chunkStart, " - ",chunkSize)
    return chunkStart, chunkSize


def DLS_EndChunk(info):
    perfInfo = [0] * 4
    tk = loc = updateJ = 0

    if info.comm == MPI.COMM_NULL or info.subChunkSize == 0:
        return
    
    tk = MPI.Wtime()

    info.t1 = tk - info.t1

    info.wStart += info.subChunkSize
    info.wSize -= info.subChunkSize
    info.sumt1 += info.t1
    info.workTime += info.t1

    if info.method == 13:
        info.sumt2 += (info.t1 * info.t1)
    
    if info.wSize == 0:
        
        if info.method == 9 or info.method == 10:
            info.mySumTimes += ((1 + info.myExecs) * info.sumt1)
            info.mySumSizes += ((1 + info.myExecs) * info.rSize)
           
        elif info.method == 11 or info.method == 12:
            info.mySumTimes += ((1 + info.myExecs) * (tk - info.t0))
            info.mySumSizes += ((1 + info.myExecs) * info.rSize)
         

        if info.method != 13:
            
            info.myIters = info.myIters + info.rSize
            info.myExecs = info.myExecs + 1
            info.sumt1 = 0.0 
            info.sumt2 = 0.0  
            info.rSize = 0

        if info.nextWRKrcvd: 
            info.t0 = MPI.Wtime()
            info.wStart = info.nextStart
            info.wSize = info.nextSize
            info.rStart = info.wStart
            info.rSize = info.wSize

            SetBreaks(info)

            info.nextSize = 0
            info.nextWRKrcvd = 0
            info.req4WRKsent = 0
    
    if info.wSize <= info.sendRequest and info.req4WRKsent==0:

        if info.method >=0 and info.method <= 8:
            perfInfo[0] = 0.0
            perfInfo[1] = 0.0
            perfInfo[2] = 1.0*info.myRank
        
        elif info.method == 9 or info.method == 10:
            nume = (info.mySumTimes + (info.myExecs+1) * info.sumt1)
            denom = (info.mySumSizes + 1.0 * (info.myExecs+1) * (info.rSize - info.wSize) )
            perfInfo[0] = nume / denom
            perfInfo[1] = 0.0
            perfInfo[2] = 1.0*info.myRank

        elif info.method == 11 or info.method == 12:
            nume = (info.mySumTimes + (info.myExecs+1) * (tk - info.t0))
            denom = (info.mySumSizes + 1.0 * (info.myExecs+1) * (info.rSize - info.wSize) )
            perfInfo[0] = nume / denom
            perfInfo[1] = 0.0
            perfInfo[2] = 1.0*info.myRank
        
        elif info.method == 13:
            perfInfo[0] = info.sumt1 / (info.rSize - info.wSize)
            
            if (info.rSize - info.wSize) > 1:
                nume = info.sumt2 - perfInfo[0] * perfInfo[0] * (info.rSize - info.wSize)
                denom = info.rSize - info.wSize - 1
                perfInfo[1] = nume / denom
                if perfInfo[1] < 0.0: 
                    perfInfo[1] = 0.0

                perfInfo[1] = math.sqrt(perfInfo[1])
            else:
                perfInfo[1] = 0.0
                perfInfo[2] = 1.0 * info.myRank

            if info.wSize == 0:
                info.myIters = info.myIters + info.rSize
                info.myExecs = info.myExecs + 1
                info.sumt1 = 0.0 
                info.sumt2 = 0.0  
                info.rSize = 0

        if info.myRank == info.foreman:

            if info.method == 9 or info.method == 10 or info.method == 11 or info.method == 12 or info.method == 13:
                loc = int(perfInfo[2])
                info.stats[3*loc+2] = info.stats[3*loc+2] + 1.0
                info.stats[3*loc] = perfInfo[0]
                info.stats[3*loc+1] = perfInfo[1]

                if info.finishedOne != info.nProcs:

                    updateJ = loc
                    for i in range(info.firstRank, info.lastRank + 1):
                        if info.stats[3 * i + 2] > 0.0 and info.stats[3 * i] < info.stats[3 * updateJ]:
                            updateJ = i
                        
                    info.finishedOne = 0

                    for i in range(info.firstRank, info.lastRank + 1):

                        if info.stats[3*i+2] == 0.0:
                            info.stats[3*i] = info.stats[3*updateJ]
                            info.stats[3*i+1] = info.stats[3*updateJ+1]
                            
                        else:
                            info.finishedOne += 1
            
            if info.chunkStart < info.lastIter:

                info.req4WRKsent = 1
                info.nextWRKrcvd = 0
                SendChunk (info, info.myRank)
            
            elif info.wSize == 0 and info.chunkStart >= info.lastIter:
                info.numENDed = info.numENDed + 1
                info.gotWork = info.numENDed != info.nProcs 

        else:
            info.comm.send(perfInfo, info.foreman, tag=REQ_TAG)
            info.req4WRKsent = 1
            info.nextWRKrcvd = 0


def  DLS_EndLoop(info):
    
    perfInfo = [0] * 4

    if info.comm == MPI.COMM_NULL:
        return

    niters = info.myIters
    worktime = info.workTime
    
    if info.method==8:
        nume = info.mySumTimes + (info.timeStep) * info.workTime
        denom = info.mySumSizes + 1.0*(info.timeStep)*(info.myIters)
        # print(nume, denom)
        perfInfo[0] = nume / denom
        perfInfo[1] = 0.0
        perfInfo[2] = 1.0*info.timeStep
        perfInfo = info.comm.gather(info.stats, root=info.foreman)

    return niters, worktime