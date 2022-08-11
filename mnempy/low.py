
import numpy as np

def transclose(phi):
    for r in range(0, phi.shape[0], 1):
        for c in range(0, phi.shape[1], 1):
            for r2 in range(0, phi.shape[0], 1):
                if phi[r][c] == 1 and phi[c][r2] == 1:
                    phi[r][r2] = 1
    np.fill_diagonal(phi, 1)
    return phi


def scorenem(D, phi, alpha, beta, rho = None, weights = None):
    if rho is scorenem.__defaults__[0]:
        rho = np.diag(np.repeat(1,D.shape[1]))
    if weights is scorenem.__defaults__[1]:
        weights = np.repeat(1, D.shape[1])
    R = D.copy()
    R[np.where(D == 1)] = np.log((1 - beta) / alpha)
    R[np.where(D == 0)] = np.log(beta / (1 - alpha))
    R = R*weights
    R = np.matmul(R,np.transpose(rho))
    phi = transclose(phi)
    score = np.matmul(R, phi)
    score = np.sum(np.max(score, axis=1))
    return score


def meanD(D, rho):
    Dm = matmul(D, np.transpose(rho))
    return Dm


def resp(D, rho, alpha, beta, phis):
    K = phis.shape[0]
    n = phis.shape[1]
    m = D.shape[0]
    p = np.reshape(np.repeat(0, K * np.shape(D)[1]), (K, np.shape(D)[1]))
    R = D.copy()
    R[np.where(D == 1)] = np.log((1 - beta) / alpha)
    R[np.where(D == 0)] = np.log(beta / (1 - alpha))
    for k in range(0, K, 1):
        theta = np.repeat(0,n*m).reshape(n,m)
        P = np.matmul(np.matmul(R,np.transpose(rho)),phis[k,:,:])
        thetaidx = np.argmax(P,1)
        for i in range(0,m,1):
            theta[thetaidx[i],i] = 1
        L = np.matmul(np.matmul(np.transpose(rho), phis[k,:,:]),np.matmul(theta,R))
        p[k,:] = np.diagonal(L)
    return p


def randphis(n, K):
    phis = np.random.randint(0, 2, K * n * n)
    phis = phis.reshape((K, n, n))
    for k in range(0, K, 1):
        phis[k, :, :] = np.triu(phis[k, :, :])
        idx = np.random.permutation(np.arange(0,n,1))
        phis[k, :, :] = phis[k, idx, :]
        phis[k, :, :] = transclose(phis[k, :, idx])
    return phis


def createdata(phis, m):
    K = phis.shape[0]
    n = phis.shape[1]
    for k in range(0, K, 1):
        if k == 0:
            D = np.repeat(np.transpose(phis[k, :, :]), m, 0)
        else:
            Dtmp = np.repeat(np.transpose(phis[k, :, :]), m, 0)
            D = np.concatenate((D, Dtmp), 1)
    return D

