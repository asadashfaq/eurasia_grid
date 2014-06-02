class FlowCalculation:
    """ This class is made for passing information about a calculation,
        to a function that solves the network structure, by only passing
        one object, thus allowing for easy multiprocessing with map from
        the multiprocessing module.
        Example: my_calc = Calculation('EU_RU_ME', 'aHE', 'q99', 'lin'),
        aHE, means that a heterogeneous array of alphas are being used,
        the optimal mixes for each region individually.

        """

    def __init__(self, layout, alphas, capacities, solvermode):
        self.layout = layout
        self.alphas = alphas
        self.capacities = capacities
        if capacities=='zerotrans':
            self.solvermode = 'raw'
        else:
            self.solvermode = solvermode


    def __str__(self):
        return ''.join([self.layout, '_', self.alphas, '_', self.capacities,
                        '_', self.solvermode])

    def copy(self):
        return FlowCalculation(self.layout, self.alphas, self.capacities,\
                                self.solvermode)

    def label(self, variationparameter):
        """ Returns a string that can be used to label a curve
            in a histgram that plots with different values of
            variationparameter.
            This is a nicer more readible string representation that str()

            """

        if variationparameter=='layout':
            if self.layout=='w':
                return 'World'
            else:
                return self.layout.replace('_','-')
        elif variationparameter=='alphas':
            if self.alphas=='aHE':
                return r'Optimal $\alpha_W$ for each region'
            if self.alphas=='aHO1':
                return r'$\alpha_W$ = 1 for all regions'
            if self.alphas=='aHO0':
                return r'$\alpha_W$ = 0 for all regions'


        elif variationparameter=='capacities':
            if self.capacities=='copper':
                return 'Unconstrained flow'
            elif self.capacities=='q99':
                return '99% quantiles'
            elif self.capacities=='hq99':
                return '0.5 * 99% quantiles'
            elif self.capacities=='zerotrans':
                return 'No transmission'
            else:
                return str(self)

        elif variationparameter=='solvermode':
            if self.solvermode=='lin':
                return 'Linear'
            elif self.solvermode=='sqr':
                return 'Square'
            elif self.solvermode=='capM':
                return 'Martin capped'
            elif self.solvermode=='capR':
                return 'Rolando capped'
            else:
                return 'No transmission'

        else:
            return str(self)

    def str_without(self, field):
        """ This function returns a string that represents the object
            but leave out one of the fields. Useful for when one parameter
            is varied, but most are kept the same, in plotting histograms
            of load and mismatch for example.

            """

        if field == 'layout':
            return ''.join([self.alphas, '_', self.capacities, '_',
                            self.solvermode])
        if field == 'alphas':
            return ''.join([self.layout, '_', self.capacities, '_',
                            self.solvermode])
        elif field == 'capacities':
            return ''.join([self.layout, '_', self.alphas, '_',
                            self.solvermode])
        elif field == 'solvermode':
            return ''.join([self.layout, '_', self.alphas, '_',
                            self.capacities])
        else:
            return str(self)
