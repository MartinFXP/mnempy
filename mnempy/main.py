
import numpy as np
from .low import transclose, scorenem, meanD, resp, randphis, createdata


def sim(K=2, n=5, m=2, s=None, alpha=0.01, beta=0.1):
    if s is sim.__defaults__[3]:
        s = K * n * 10
    phis = randphis(n, K)
    sidx = np.random.randint(0, K * n, s)
    rho = np.repeat(0, n * sidx.__len__()).reshape(n, sidx.__len__())
    for i in range(0, sidx.__len__()):
        rho[sidx[i] % n, i] = 1
    D = createdata(phis, m)
    D = D[:, sidx]
    ones = np.where(D == 1)
    zeros = np.where(D == 0)
    alpha = 0.01
    beta = 0.1
    fp = np.random.choice(zeros[0].__len__(), int(np.floor(alpha * zeros[0].__len__()) + 1))
    fn = np.random.choice(ones[0].__len__(), int(np.floor(beta * ones[0].__len__()) + 1))
    D[zeros[0][fp], zeros[1][fp]] = 1
    D[ones[0][fn], ones[1][fn]] = 0
    output = {
        "D": D,
        "phis": phis,
        "rho": rho
    }
    return output


def nem(D, alpha = 0.01, beta = 0.1, step = 0.01, rho = None, weights = None):
    if rho is nem.__defaults__[3]:
        rho = np.diag(np.repeat(1,D.shape[1]))
    if weights is nem.__defaults__[4]:
        weights = np.repeat(1, D.shape[1])
    n = rho.shape[0]
    phi = np.zeros((n, n))
    phi = transclose(phi)
    ll = scorenem(D, phi, alpha, beta, rho, weights)
    llold = ll - 1
    while ll > llold:
        llold = ll
        ibest = n + 1
        jbest = n + 1
        besttype = 0
        llbest = ll
        for i in range(0, n, 1):
            for j in range(0, n, 1):
                phinew = phi.copy()
                phinew[i][j] = 1 - phinew[i][j]
                ll = scorenem(D, phinew, alpha, beta, rho, weights)
                if ll > llbest:
                    llbest = ll
                    ibest = i
                    jbest = j
        alphanew = min(alpha + step, 1 - step)
        ll = scorenem(D, phi, alphanew, beta, rho, weights)
        if ll > llbest:
            llbest = ll
            besttype = 1
        alphanew = max(alpha - step, step)
        ll = scorenem(D, phi, alphanew, beta, rho, weights)
        if ll > llbest:
            llbest = ll
            besttype = 2
        betanew = min(beta + step, 1 - step)
        ll = scorenem(D, phi, alpha, betanew, rho, weights)
        if ll > llbest:
            llbest = ll
            besttype = 3
        betanew = max(beta - step, step)
        ll = scorenem(D, phi, alpha, betanew, rho, weights)
        if ll > llbest:
            llbest = ll
            besttype = 4
        if ibest <= n and jbest <= n and besttype == 0:
            phi[ibest][jbest] = 1 - phi[ibest][jbest]
        if besttype == 1:
            alpha = alpha + step
        if besttype == 2:
            alpha = alpha - step
        if besttype == 3:
            beta = beta + step
        if besttype == 4:
            beta = beta - step
        ll = llbest
    output = {
        "phi": phi,
        "lh": ll,
        "alpha": alpha,
        "beta": beta
    }
    return output


def mnem(D, rho, K = 2, alpha = 0.01, beta = 0.1, step = 0.01):
    n = rho.shape[0]
    l = D.shape[1]
    pi = np.random.uniform(0.1,0.9,K)
    pi = pi/pi.sum()
    phis = randphis(n, K)
    p = resp(D,rho,alpha,beta,phis)
    lh = np.sum(np.log(np.matmul(np.repeat(1,K).reshape(1,K),np.exp(p)*np.repeat(pi,l).reshape(K,l))))
    lhold = lh - 1
    while lh > lhold:
        probs = np.transpose(np.transpose(np.exp(p))*pi)
        probs = probs/np.sum(probs,0)
        lhold = lh
        for k in range(0,K,1):
            result = nem(D,alpha,beta,step,rho,weights=probs[k,:])
            phis[k,:,:] = result["phi"]
        p = resp(D, rho, alpha, beta, phis)
        lh = np.sum(np.log(np.matmul(np.repeat(1,K).reshape(1,K),np.exp(p)*np.repeat(pi,l).reshape(K,l))))
    output = {
        "phis": phis,
        "lh": lh,
        "resp": p
    }
    return output

