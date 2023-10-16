import conditional_independence
import numpy as np
import pandas as pd
import scipy.linalg as slin
import scipy.optimize as sopt
import sys

from causaldag import igsp

sys.path.append("../python_methods_comparison/golem/src/")
sys.path.append("../python_methods_comparison/golem/src/models")
sys.path.append("../python_methods_comparison/golem/src/trainers")

import golem

def fit_golem(X, lam):
    # Do I need to copy?
    W0 = golem.golem(X, lambda_1=lam, lambda_2=5.0,
                     equal_variances=True, seed=np.random.randint(0, 2**32-1))
    W = golem.golem(X, lambda_1=lam, lambda_2=5.0,
                    equal_variances=False, seed=np.random.randint(0, 2**32-1), B_init=W0)
    return W


def loss_notears(W, X, loss_type='l2'):
    """Evaluate value and gradient of loss."""
    M = X @ W
    if loss_type == 'l2':
        R = X - M
        loss = 0.5 / X.shape[0] * (R ** 2).sum()
        G_loss = - 1.0 / X.shape[0] * X.T @ R
    elif loss_type == 'logistic':
        loss = 1.0 / X.shape[0] * (np.logaddexp(0, M) - X * M).sum()
        G_loss = 1.0 / X.shape[0] * X.T @ (sigmoid(M) - X)
    elif loss_type == 'poisson':
        S = np.exp(M)
        loss = 1.0 / X.shape[0] * (S - X * M).sum()
        G_loss = 1.0 / X.shape[0] * X.T @ (S - X)
    else:
        raise ValueError('unknown loss type')
#         print(loss, G_loss)
    return loss, G_loss	


def notears_linear(X, lambda1, loss_type='l2', max_iter=100, h_tol=1e-8, rho_max=1e+16, w_threshold=0.3):
    """Solve min_W L(W; X) + lambda1 ‖W‖_1 s.t. h(W) = 0 using augmented Lagrangian.
    Args:
        X (np.ndarray): [n, d] sample matrix
        lambda1 (float): l1 penalty parameter
        loss_type (str): l2, logistic, poisson
        max_iter (int): max num of dual ascent steps
        h_tol (float): exit if |h(w_est)| <= htol
        rho_max (float): exit if rho >= rho_max
        w_threshold (float): drop edge if |weight| < threshold
    Returns:
        W_est (np.ndarray): [d, d] estimated DAG
    """
    def _loss(W):
        """Evaluate value and gradient of loss."""
        M = X @ W
        if loss_type == 'l2':
            R = X - M
            loss = 0.5 / X.shape[0] * (R ** 2).sum()
            G_loss = - 1.0 / X.shape[0] * X.T @ R
        elif loss_type == 'logistic':
            loss = 1.0 / X.shape[0] * (np.logaddexp(0, M) - X * M).sum()
            G_loss = 1.0 / X.shape[0] * X.T @ (sigmoid(M) - X)
        elif loss_type == 'poisson':
            S = np.exp(M)
            loss = 1.0 / X.shape[0] * (S - X * M).sum()
            G_loss = 1.0 / X.shape[0] * X.T @ (S - X)
        else:
            raise ValueError('unknown loss type')
            
#         print(loss, G_loss)
        return loss, G_loss

    def _h(W):
        """Evaluate value and gradient of acyclicity constraint."""
        E = slin.expm(W * W)  # (Zheng et al. 2018)
        h = np.trace(E) - d
        #     # A different formulation, slightly faster at the cost of numerical stability
        #     M = np.eye(d) + W * W / d  # (Yu et al. 2019)
        #     E = np.linalg.matrix_power(M, d - 1)
        #     h = (E.T * M).sum() - d
        G_h = E.T * W * 2
        return h, G_h

    def _adj(w):
        """Convert doubled variables ([2 d^2] array) back to original variables ([d, d] matrix)."""
        return (w[:d * d] - w[d * d:]).reshape([d, d])

    def _func(w):
        """Evaluate value and gradient of augmented Lagrangian for doubled variables ([2 d^2] array)."""
        W = _adj(w)
        loss, G_loss = _loss(W)
        h, G_h = _h(W)
        obj = loss + 0.5 * rho * h * h + alpha * h + lambda1 * w.sum()
        G_smooth = G_loss + (rho * h + alpha) * G_h
        g_obj = np.concatenate((G_smooth + lambda1, - G_smooth + lambda1), axis=None)
        return obj, g_obj

    n, d = X.shape
    w_est, rho, alpha, h = np.zeros(2 * d * d), 1.0, 0.0, np.inf  # double w_est into (w_pos, w_neg)
    bnds = [(0, 0) if i == j else (0, None) for _ in range(2) for i in range(d) for j in range(d)]
    if loss_type == 'l2':
        X = X - np.mean(X, axis=0, keepdims=True)
    for i in range(max_iter):
        w_new, h_new = None, None
        while rho < rho_max:
            sol = sopt.minimize(_func, w_est, method='L-BFGS-B', jac=True, bounds=bnds)
            w_new = sol.x
            h_new, _ = _h(_adj(w_new))
            if h_new > 0.25 * h:
                rho *= 10
            else:
                break
        w_est, h = w_new, h_new
        alpha += rho * h
        # print(i, sol.fun, h)
        if h <= h_tol or rho >= rho_max:
            break
    W_est = _adj(w_est)
    W_est[np.abs(W_est) < w_threshold] = 0
    return W_est


def scale_data(data):
    scaler = StandardScaler()
    data_scaled = {}
    for k, v in data.items():
        scaled = scaler.fit_transform(v)
        data_scaled[k] = scaled
        
    return data_scaled


class DOTEARS:
    def __init__(self, data, lambda1=0.1, loss_type='l2', max_iter=100, h_tol=1e-8, rho_max=1e+16, scaled=False, w_threshold=0, obs_only=False):
        self.data = data
        self.lambda1 = lambda1
        self.loss_type = loss_type
        self.max_iter = max_iter
        self.h_tol = h_tol
        self.rho_max = rho_max
        self.w_threshold = w_threshold
        self.scaled = scaled
        self.obs_only = obs_only

        self.p = data['obs'].shape[1]
        self.V_inverse = (self.estimate_exogenous_variances(data) ** (-1)) * np.identity(self.p)
        
        self.scaled = scaled
        if self.scaled:
#             self.original_variances = {}
#             for k, v in data.items():
#                 self.original_variances[k] = v.var(axis=0)
            self.data = scale_data(data)
    
    def estimate_exogenous_variances(self, data):
        p = data['obs'].shape[1]
        variances = np.zeros(p)
        for k, v in data.items():
            if k == 'obs':
                continue
            variances[int(k)] = v.var(axis=0)[int(k)]

        return variances

    def loss(self, W):
        data = self.data
        p = self.p
        obs_only = self.obs_only

        V_inverse = self.V_inverse
        V_half = V_inverse ** 0.5

        if self.loss_type == 'l2':
            G_loss = 0
            loss = 0

            for j in data.keys():
                if self.obs_only:
                    if j != 'obs':
                        continue

                mask = np.ones((p, p))

                if j != 'obs':
                    mask[:, int(j)] = 0

                W_j = mask * W

                R = data[j] - data[j] @ W_j
#                 if self.scaled:
 #                    V_half = (V_inverse * self.original_variances[j]) ** 0.5

                # new gradient
                loss += 0.5 / data[j].shape[0] * ((R @ V_half) ** 2).sum()
#                 if self.scaled:
#                     G_loss += - 1.0 / data[j].shape[0] * data[j].T @ R @ (V_inverse * self.original_variances[j])
#                 else:
                G_loss += - 1.0 / data[j].shape[0] * data[j].T @ R @ V_inverse

            if not self.obs_only:
                loss /= len(data.keys())
                G_loss /= len(data.keys())
        else:
            raise ValueError('unknown loss type')
        return loss, G_loss

    def _h(self, W):
        """Evaluate value and gradient of acyclicity constraint."""
        E = slin.expm(W * W)  # (Zheng et al. 2018)
        h = np.trace(E) - self.p
        #     # A different formulation, slightly faster at the cost of numerical stability
        #     M = np.eye(d) + W * W / d  # (Yu et al. 2019)
        #     E = np.linalg.matrix_power(M, d - 1)
        #     h = (E.T * M).sum() - d
        G_h = E.T * W * 2
        return h, G_h

    def _adj(self, w):
        """Convert doubled variables ([2 d^2] array) back to original variables ([d, d] matrix)."""
        d = self.p
        return (w[:d * d] - w[d * d:]).reshape([d, d])

    def _func(self, w):
        """Evaluate value and gradient of augmented Lagrangian for doubled variables ([2 d^2] array)."""
        rho = self.rho
        W = self._adj(w)
#         W = np.triu(W, 1) # added
        loss, G_loss = self.loss(W)
        h, G_h = self._h(W)
        obj = loss + 0.5 * rho * h * h + self.alpha * h + self.lambda1 * w.sum()
        G_smooth = G_loss + (rho * h + self.alpha) * G_h
        g_obj = np.concatenate((G_smooth + self.lambda1, - G_smooth + self.lambda1), axis=None)
        return obj, g_obj

    def fit(self):
        d = self.p
        data = self.data
        w_est, self.rho, self.alpha, h = np.zeros(2 * d * d), 1.0, 0.0, np.inf  # double w_est into (w_pos, w_neg)

        bnds = [(0, 0) if i == j else (0, None) for _ in range(2) for i in range(d) for j in range(d)]
        for k, v in data.items():
            data[k] = data[k] - np.mean(data[k], axis=0, keepdims=True)

        for iter_number in range(self.max_iter):
            w_new, h_new = None, None
            while self.rho < self.rho_max:
                sol = sopt.minimize(self._func, w_est, method='L-BFGS-B', jac=True, bounds=bnds)
                w_new = sol.x
                h_new, _ = self._h(self._adj(w_new))

                if h_new > 0.25 * h:
                    self.rho *= 10
                else:
                    break    
            w_est, h = w_new, h_new
            self.alpha += self.rho * h
            # print(iter_number, sol.fun, h)
            if h <= self.h_tol or self.rho >= self.rho_max:
                break
        W_est = self._adj(w_est)
        W_est[np.abs(W_est) < self.w_threshold] = 0
        return W_est


def make_dotears_data(X, targets):
  data_dotears = {}
  for key in np.unique(targets):
    if key != "control":
      data_dotears[int(key.replace("V", ""))-1] = X[targets==key,]
    else:
      data_dotears["obs"] = X[targets==key,]
  return data_dotears


# Adapted from https://github.com/asxue/dotears/blob/main/workflow/scripts/igsp.py
def run_igsp(X_dotears, alpha=0.001, alpha_inv=0.001):
    obs_data = X_dotears['obs']
    p = obs_data.shape[1]
    nodes = list(range(p))

    # iv_samples_list is a list of n x p ndarrays
    inv_samples_from_data = [v for k, v in X_dotears.items() if k != 'obs']

    # setting_list is list of dicts
    # each dict is key 'intervention' to a list of nodes
    settings_from_data = [dict(interventions=[k]) for k, v in X_dotears.items() if k != 'obs']

    obs_suffstat = conditional_independence.partial_correlation_suffstat(obs_data)
    invariance_suffstat = conditional_independence.gauss_invariance_suffstat(
        obs_data, inv_samples_from_data)

    ci_tester = conditional_independence.MemoizedCI_Tester(
        conditional_independence.partial_correlation_test, obs_suffstat, alpha=alpha)
    invariance_tester = conditional_independence.MemoizedInvarianceTester(
        conditional_independence.gauss_invariance_test, invariance_suffstat, alpha=alpha_inv)

    est_dag = igsp(settings_from_data, nodes, ci_tester, invariance_tester)
    W_igsp = est_dag.to_amat()[0]
    return(W_igsp)
