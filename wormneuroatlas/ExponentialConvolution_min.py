import numpy as np

class ExponentialConvolution_min:
    
    def __init__(self,exp):
        self.exp = [exp]
        
    def __call__(self,time):
        return self.eval(time)
            
    def eval(self,x,dtype=np.float64):
        '''Evaluates the ExponentialConvolution in the time domain.
        
        Parameters
        ----------
        x: array_like
            Time axis. All times should be positive.
        dtype: type (optional)
            Type of the output array. Default: np.float64        
        drop_branches: int or array_like of int
            Branches to be ignored in the evaluation. Default: None.
            
        Returns
        -------
        out: numpy.ndarray
            ExponentialConvolution evaluated on x.
        '''
        assert np.all(x>=0)
        
        
        out = np.zeros_like(x,dtype=dtype)
            
        for exp in self.exp[-1]:
            # Skip terms that are in excluded branches
            branch = exp["branch"]
                
            g = exp["g"]
            factor = exp["factor"]
            power_t = exp["power_t"]
            
            if power_t==0: mult=1.
            else: mult=np.power(x,power_t)
            out += factor*mult*np.exp(-g*x)
            
        return out
    
    @staticmethod
    def eval2(x,gs,cs,pt,dtype=np.float64):
        '''Static alternative evaluation method that takes the parameters of
        ExponentialConvolution as inputs (obtainable from 
        ExponentialConvolution.get_bare_params()).
        
        Parameters
        ----------
        x: array_like
            Time axis. All times should be positive.
        gs: array_like of floats
            Gammas
        cs: array_like of floats
            Amplitudes
        pt: array_like of floats
            Powers of t
            
        Returns
        -------
        out: numpy.ndarray
            ExponentialConvolution evaluated on x.
        '''
        assert np.all(x>=0)
        out = np.zeros_like(x,dtype=dtype)
            
        for i in np.arange(len(gs)):
            g = gs[i]
            factor = cs[i]
            power_t = pt[i]
            if power_t==0: mult=1.
            else: mult=np.power(x,power_t)
            out += factor*mult*np.exp(-g*x)
            
        return out
        
    def multiply_scalar_inplace(self,factor):
        for i in np.arange(len(self.exp[-1])):
            self.exp[-1][i]["factor"]*=factor
