#! /usr/bin/env python3

import numpy as np
import lmfit
import klassez as kz

#fitting functions from TRAGICO

def exponential_ls(param, x, y, multi=1, result=False):
    """
    Least squares residuals for an exponential function with multiplicity from 1 to 3, 
    with the option to return the model.
    
    .. math::
    
        model = A*f(x, k) + a, 
    
    
    where :math:`f(x, k)` is the exponential function with multiplicity `multi` and parameters `k`, :math:`A` is the optimal scaling factor and :math:`a` is the optimal offset in the least squares sense.
    The model is computed by :func:`exponential_model`, which is called by this function.
    This function is taken from the `TRAGICO code`_.
    
    .. _TRAGICO code: https://github.com/letiziafiorucci/tragico  
    
    Parameters
    ----------
    param : lmfit.Parameters
        The parameters for the exponential model.
    x : array_like
        The independent variable data.
    y : array_like
        The dependent variable data.
    multi : int, optional
        The multiplicity of the exponential model (1, 2, or 3). Default is 1.
    result : bool, optional
        If True, return the model and the optimal A and a in the LS sense. Default is False.

    Returns
    -------
    array_like
        The residuals of the exponential model.
    or
    tuple
        If result is True, returns a tuple containing the model, optimal A, and optimal a.
    """
    
    model = exponential_model(param, x, multi)
    #print(model)
    den = np.mean(model**2)-np.mean(model)**2
    a = (np.mean(model**2)*np.mean(y)-np.mean(model*y)*np.mean(model))/den
    A = (np.mean(model*y)-(np.mean(model)*np.mean(y)))/den
  
    model *= A
    model += a  

    if any(np.isnan(model)):
        print(f'Nan in {multi}exp model')
    if result==False:

        return y-model
    else:
        return model, A, a
    
def exponential_model(param, x, multi=1, A=1, a=0):
    """
    Exponential model with multiplicity from 1 to 3.
    
    .. math::
    
        model = A*f(x, k) + a, 
    
    
    where :math:`f(x, k)` is the exponential function with multiplicity `multi` and parameters `k`, :math:`A` is the optimal scaling factor and :math:`a` is the optimal offset in the least squares sense.
    This function is taken from the `TRAGICO code`_.
    
    .. _TRAGICO code: https://github.com/letiziafiorucci/tragico  
    
    Parameters
    ----------
    param : lmfit.Parameters
        The parameters for the exponential model.
    x : array_like
        The independent variable data.
    multi : int, optional
        The multiplicity of the exponential model (1, 2, or 3). Default is 1.
    A : float, optional
        The scaling factor. Default is 1.
    a : float, optional
        The offset. Default is 0.

    Returns
    -------
    array_like
        The exponential model.
    """
    par=param.valuesdict()
    if multi==1:
        k = np.exp(par['k'])
        model = np.exp(-x*k)
    elif multi==2:
        k1 = np.exp(par['k1'])
        k2 = np.exp(par['k2'])
        f1 = 1 / (1 + np.exp(-par['f1']))  # Sigmoid transformation: maps (-∞, ∞) to (0, 1)
        model = f1*np.exp(-x*k1)+(1-f1)*np.exp(-x*k2)
    elif multi==3:
        k1 = np.exp(par['k1'])
        k2 = np.exp(par['k2'])
        k3 = np.exp(par['k3'])
        f1 = 1 / (1 + np.exp(-par['f1']))  # Sigmoid transformation: maps (-∞, ∞) to (0, 1)
        f2 = 1 / (1 + np.exp(-par['f2']))  # Sigmoid transformation: maps (-∞, ∞) to (0, 1)
        f2 = f2 * (1 - f1)                  # Maps (0, 1) to (0, 1-f1)
        model = f1*np.exp(-x*k1)+f2*np.exp(-x*k2)+(1-f1-f2)*np.exp(-x*k3)

    model *= A
    model += a

    return model
    
def fit_exponential(x, y, multi=1):
    """
    Fit an exponential function with multiplicity from 1 to 3 to the data (x, y) using least squares optimization.
    Calls :func:`exponential_ls` to compute the residuals and :func:`exponential_model` to compute the model.

    Parameters
    ----------
    x : array_like
        The independent variable data.
    y : array_like
        The dependent variable data.
    multi : int, optional
        The multiplicity of the exponential model (1, 2, or 3). Default is 1.

    Returns
    -------
    lmfit.MinimizerResult
        The result of the least squares optimization.
    """
    param = lmfit.Parameters()
    if multi==1:
        param.add('k', value=0, min=-5, max=5)
    if multi==2:
        param.add('k1', value=0, min=-5, max=5)
        param.add('k2', value=0, min=-5, max=5)
        param.add('f1', value=0, min=-5, max=5)  # f1 will be transformed to (0, 1) in the fitting function
    if multi==3:
        param.add('k1', value=0, min=-5, max=5)
        param.add('k2', value=0, min=-5, max=5)
        param.add('k3', value=0, min=-5, max=5)
        param.add('f1', value=0, min=-5, max=5)  # f1 will be transformed to (0, 1) in the fitting function
        param.add('f2', value=0, min=-5, max=5)  # f2 will be transformed to (0, 1-f1) in the fitting function
    minner = lmfit.Minimizer(exponential_ls, param, fcn_args=(x, y), fcn_kws={'multi': multi})
    result = minner.minimize()
    return result

#Functions for fitting spectra envelopes

def fit_skewnormal(x,y):
    r"""
    Fits the NH region of the 1D spectrum to a skew normal distribution using least squares optimization. The function returns the result of the optimization, which contains the fitted parameters of the skew normal distribution. The skew normal distribution is defined as:
    
    .. math:: 
    
    f(x) = A \cdot \frac{1}{\sigma \sqrt{2\pi}} e^{-\frac{(x-\mu)^2}{2\sigma^2}} \left(1 + \text{erf}\left(\alpha \frac{x-\mu}{\sigma \sqrt{2}}\right)\right) + a
    
    Parameters
    ----------
    x : array_like
        The independent variable data (ppm values).
    y : array_like
        The dependent variable data (intensity values).
        
    Returns
    -------
    lmfit.MinimizerResult
        The result of the least squares optimization, which contains the fitted parameters of the skew normal distribution
    """
    param = lmfit.Parameters()
    param.add('a', value=4)
    param.add('u', value=8.25)
    param.add('s', value=0.6)
    minner = lmfit.Minimizer(skgaussian_ls, param, fcn_args=(x, y))
    result = minner.minimize()
    return result

def skgaussian_ls(param, x, y, result=False):
    """
    Computes the skew Gaussian model and optionally returns the model along with the scaling parameters.

    Parameters
    ----------
    param : lmfit.Parameters
        The parameters for the skew Gaussian model.
    x : array_like
        The independent variable data (ppm values).
    y : array_like
        The dependent variable data (intensity values).
    result : bool, optional
        If True, returns the model along with the scaling parameters. Default is False.

    Returns
    -------
    array_like
        The residuals (y - model) if result is False.
    tuple
        The model, A, and a if result is True.
    """
    
    model = kz.sim.f_skgaussian(x, param['u'].value, param['s'].value, 1, param['a'].value)
    #print(model)
    den = np.mean(model**2)-np.mean(model)**2
    a = (np.mean(model**2)*np.mean(y)-np.mean(model*y)*np.mean(model))/den
    A = (np.mean(model*y)-(np.mean(model)*np.mean(y)))/den
  
    model *= A
    model += a  

    if any(np.isnan(model)):
        print(f'Nan in skgaussian model')
    if result==False:

        return y-model
    else:
        return model, A, a
