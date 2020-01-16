######################################
#   Basic normalization functions
######################################

def cpm(X):
    X_norm = np.array([
        x/sum(x)
        for x in X
    ])
    X_norm *= 10e6
    return X_norm


def log_cpm(X):
    X_cpm = cpm(X)
    return np.log(X_cbpm+1)
