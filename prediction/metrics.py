"""

"""

from numpy import mean, median, abs, sum, cumsum, histogram, sqrt
from sklearn.metrics import mean_squared_error, r2_score



def rmse(y_true, y_pred):
    """
rmse
    description:
        computes RMSE
"""
    return sqrt(mean_squared_error(y_true, y_pred))


def all_metrics(y_true, y_pred):
    """
metrics
    description:
        computes a standard set of performance metrics using true and predicted values
            * training and test set R-squared (R2)
            * training and test set root mean squared error (RMSE)
            * training and test set mean absolute error (MAE)
            * cumulative error distribution at the <1, <3, <5, and <10% levels
              for training and test set (CE135A)
"""   
    # compute metrics
    abs_y_err = abs(y_pred - y_true)
    r2 = r2_score(y_true, y_pred)
    mae = mean(abs_y_err)
    mdae = median(abs_y_err)
    rmse_ = rmse(y_true, y_pred)
    y_err_percent = 100. * abs_y_err / y_true
    mre = mean(y_err_percent)
    mdre = median(y_err_percent)
    cum_err = cumsum(histogram(y_err_percent, [_ for _ in range(101)])[0])
    cum_err = 100. * cum_err / sum(cum_err)
    ce1, ce3, ce5, ceA = cum_err[0], cum_err[2], cum_err[4], cum_err[9]
    summary = {'R2': r2, 'MAE': mae, 'MDAE': mdae, 'MRE': mre, 'MDRE': mdre, 'RMSE': rmse_, 'CE135A': [ce1, ce3, ce5, ceA]}
    return summary


