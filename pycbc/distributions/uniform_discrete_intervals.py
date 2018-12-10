import numpy
from pycbc.distributions import uniform
from pycbc.distributions import bounded

class UniformIntervals(bounded.BoundedDist):
    name = 'uniform_intervals'

    def __init__(self, **params):
         # temporarily suppress numpy divide by 0 warning
        self._stride = {}
        self._width = {}

        stride_args = [p for p in params if p.endswith('_stride')]
        width_args = [p for p in params if p.endswith('_width')]

        self._stride = dict([[p.split("_stride")[0], params.pop(p)] for p in stride_args])
        self._width = dict([[p.split("_width")[0], params.pop(p)] for p in width_args])

        uni_distr_obj = super(UniformIntervals, self).__init__(**params)
        missing = set(self._stride.keys(), self._width.keys()) - set(params.keys())

        if any(missing):
            raise ValueError("stride provided for unknown params {}".format(
                             ', '.join(missing)))
        self._stride.update(dict([[p, 0.]
            for p in params if p not in self._stride]))

        self._width.update(dict([[p, 0.]
            for p in params if p not in self._width]))

    @property
    def norm(self, size=1):
        return size * self._norm

    @property
    def lognorm(self, size=1):
        return numpy.log(norm(size=size))

    @property
    def stride(self):
        return self._stride

    def _pdf(self, size=1, param=None, **kwargs):
        """Returns the pdf at the given values. The keyword arguments must
        contain all of parameters in self's params. Unrecognized arguments are
        ignored.
        """
        err_str = "The pdf of this distribution is not well defined! Returning 1.0."
        raise NotImplementedError(err_str)
        return 1.0

    def _logpdf(self, size=1, **kwargs):
        """Returns the log of the pdf at the given values. The keyword
        arguments must contain all of parameters in self's params. Unrecognized
        arguments are ignored.
        """
        err_str = "The log pdf of this distribution is not well defined! Returning 0.0."
        raise NotImplementedError(err_str)
        return 0.0

    def rvs(self, size=1, param=None):
        """Gives a set of random values drawn from this distribution.
        Parameters
        ----------
        size : {1, int}
            The number of values to generate; default is 1.
        param : {None, string}
            If provided, will just return values for the given parameter.
            Otherwise, returns random values for each parameter.
        Returns
        -------
        structured array
            The random values in a numpy structured array. If a param was
            specified, the array will only have an element corresponding to the
            given parameter. Otherwise, the array will have an element for each
            parameter in self's params.
        """
        if param is not None:
            dtype = [(param, float)]
        else:
            dtype = [(p, float) for p in self.params]

        arr = numpy.zeros(size, dtype=dtype)
        for (p,_) in dtype:
            x = numpy.arange(0, size, 1)
            b = self.bounds[p][0] + self.stride[p] * x + self.width[p]
            a = self.bounds[p][0] + self.stride[p] * x
            arr[p] = numpy.random.uniform(a, b) 

        for (p_) in dtype:
            if arr[p] > self.bounds[p][1]:
                raise UserWarning("Caution, some samples lie beyond the bounds. Select different " \
                                  + "stride, width, number of samples, or bounds.")

        return arr

    @classmethod
    def from_config(cls, cp, section, variable_args):
        """Returns a distribution based on a configuration file. The parameters
        for the distribution are retrieved from the section titled
        "[`section`-`variable_args`]" in the config file.
        Parameters
        ----------
        cp : pycbc.workflow.WorkflowConfigParser
            A parsed configuration file that contains the distribution
            options.
        section : str
            Name of the section in the configuration file.
        variable_args : str
            The names of the parameters for this distribution, separated by
            ``VARARGS_DELIM``. These must appear in the "tag" part
            of the section header.
        Returns
        -------
        Uniform
            A distribution instance from the pycbc.inference.prior module.
        """
        return super(UniformIntervals, cls).from_config(cp, section, variable_args,
                     bounds_required=True)


__all__ = ['UniformIntervals']
