import math
import sys

OPT_TOL = 0.0001  # tolerance for iterative optimization


def EMUpdate(meansIn, te_features, uniq_counts, multi_reads, estimatedReadLength):
    # reassign multi-reads proportionally to the relative abundance of each TE
    meansOut = [0] * len(meansIn)

    multi_counts = computeAbundances(meansIn, multi_reads)
    sys.stderr.write("total multi counts = " + str(sum(multi_counts)) + "\n")
    for tid in range(len(meansIn)):
        tlen = te_features.getLength(tid)
        if tlen < 0:
            sys.stderr.write("Error in optimization: the TE does not exist!\n")
            raise
        effectiveLength = tlen - estimatedReadLength + 1
        if effectiveLength > 0:
            meansOut[tid] = (uniq_counts[tid] + multi_counts[tid]) / effectiveLength  # py3
        else:
            meansOut[tid] = 0.0

    summeans = sum(meansOut)
    meansOut = [float(i)/summeans for i in meansOut]

    return meansOut


def EMestimate(te_features, multi_reads, uniq_counts, multi_counts, numItr, estimatedReadLength):
    sys.stderr.write("multi-reads = %s " % str(len(multi_reads)))
    means0 = []
    for tid in range(len(uniq_counts)):
        tlen = te_features.getLength(tid)
        if tlen < 0:
            sys.stderr.write("Error in optimization: the TE does not exist!\n")
            raise
        effectiveLength = tlen - estimatedReadLength + 1
        if effectiveLength > 0:
            means0.append(1.0 * (uniq_counts[tid] + multi_counts[tid]) / effectiveLength)
        else:
            means0.append(0.0)

    # relative abundance
    summeans0 = sum(means0)
    means0 = [float(i)/summeans0 for i in means0]
    '''
    /**
         * Defaults for these values taken from the R implementation of
         * [SQUAREM](http://cran.r-project.org/web/packages/SQUAREM/index.html).
         */
    '''
    minStep = 1.0
    maxStep = 1.0
    mStep = 4.0
    cur_iter = 0
    t_size = len(uniq_counts)
    r = [0] * t_size
    r2 = [0] * t_size
    v = [0] * t_size
    meansPrime = [0.0] * t_size
    while cur_iter < numItr:
        cur_iter += 1
        sys.stderr.write("SQUAREM iteraton [" + str(cur_iter) + "]\n")
        try:
            means1 = EMUpdate(means0, te_features, uniq_counts, multi_reads, estimatedReadLength)
        except:
            sys.stderr.write("Error in EMupdate\n")
            raise

        try:
            means2 = EMUpdate(means1, te_features, uniq_counts, multi_reads, estimatedReadLength)
        except:
            sys.stderr.write("Error in EMupdate\n")
            raise

        for tid in range(len(means0)):
            r[tid] = means1[tid] - means0[tid]
            r2[tid] = means2[tid] - means1[tid]
            v[tid] = (means2[tid] - means1[tid]) - r[tid]

        rNorm = math.sqrt(sum([a * b for a, b in zip(r, r)]))
        r2Norm = math.sqrt(sum([a * b for a, b in zip(r2, r2)]))
        vNorm = math.sqrt(sum([a * b for a, b in zip(v, v)]))
        rr = sum([a * b for a, b in zip(r, v)])
        rvNorm = math.sqrt(abs(rr))

        if vNorm == 0:
            means0 = means1
            sys.stderr.write("at iteration " + str(cur_iter) + " vNorm == 0 \n")
            break
        alphaS = rNorm / rvNorm
        alphaS = max(minStep, min(maxStep, alphaS))

        if rNorm < OPT_TOL:
            sys.stderr.write("rNome < OPT_TOL \n")
            break
        if r2Norm < OPT_TOL:
            sys.stderr.write("r2Nome < OPT_TOL \n")
            means0 = means2
            break

        for tid in range(len(means0)):
            meansPrime[tid] = max(0.0, means0[tid] + 2 * alphaS * r[tid] + (alphaS * alphaS) * v[tid])

        # Stabilization step
        if abs(alphaS - 1.0) > 0.01:
            sys.stderr.write("alpha = " + str(alphaS) + ".\n ")
            sys.stderr.write("Performing a stabilization step.\n")
            try:
                meansPrime = EMUpdate(meansPrime, te_features, uniq_counts, multi_reads, estimatedReadLength)
            except:
                sys.stderr.write("Error in EMupdate\n")
                raise

        sys.stderr.write("alpha = " + str(alphaS) + ", ")

        if alphaS == maxStep:
            maxStep = mStep * maxStep

        if 0 > minStep == alphaS:
            minStep = mStep * minStep

        meansPrime, means0 = means0, meansPrime

    if cur_iter >= numItr:
        sys.stderr.write("did not converge by " + str(numItr) + " iterations. \n")
    else:
        sys.stderr.write("converge at iteration " + str(cur_iter) + "\n")
    new_multi_counts = computeAbundances(means0, multi_reads)
    return new_multi_counts


def computeAbundances(meansIn, multi_reads):
    size = len(meansIn)
    multi_counts = [0] * size

    sys.stderr.write("num of multi reads = " + str(len(multi_reads)) + "\n")
    for kid in range(len(multi_reads)):
        TE_transcripts = multi_reads[kid]
        totalMass = 0.0
        for tid in TE_transcripts:
            totalMass += meansIn[tid]

        if totalMass > 0.0:
            norm = 1.0 / totalMass
        else:
            norm = 0.0

        for tid in TE_transcripts:
            multi_counts[tid] += meansIn[tid] * norm

    return multi_counts
