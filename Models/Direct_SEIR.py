import matplotlib.pyplot as plt

class SEIRModel:
    
    def __init__(self, exposure_rate=1.0, incubation_rate=1.0, recovery_rate=0.1, sus0 = 0.99, exp0 = 0.01, inf0=0.0, rec0=0.0, days=100): ##initialised parameters
        self.exposure_rate = exposure_rate
        self.incubation_rate = incubation_rate
        self.recovery_rate = recovery_rate
        self.days = days
        self.__days_list = range(0,days + 1) ##list of days
        self.__sus_list = [sus0] 
        self.__exp_list = [exp0]
        self.__inf_list = [inf0]
        self.__rec_list = [rec0] ##initial list of population

    def __update_sus (self, sus, inf):
        return sus - self.exposure_rate * inf * sus

    def __update_exp (self, exp, sus, inf):
        return exp + self.exposure_rate * inf * sus - self.incubation_rate * exp 

    def __update_inf (self, inf, exp):
        return inf + self.incubation_rate * exp - self.recovery_rate * inf

    def __update_rec (self, rec, inf):
        return rec + self.recovery_rate * inf ##calculate the rates and returns updated population size

    def calculate(self):
        for i in range(0, self.days):
            sus = self.__sus_list[-1]
            exp = self.__exp_list[-1]
            inf = self.__inf_list[-1]
            rec = self.__rec_list[-1]
            
            self.__sus_list.append(self.__update_sus(sus, inf))
            self.__exp_list.append(self.__update_exp(exp, sus, inf))
            self.__inf_list.append(self.__update_inf(inf, exp))
            self.__rec_list.append(self.__update_rec(rec, inf)) ##updates population lists of each agent

    def plot(self):
        fig, ax = plt.subplots()
        plt.plot(self.__days_list, self.__sus_list, label = "Susceptible")
        plt.plot(self.__days_list, self.__exp_list, label = "Exposed")
        plt.plot(self.__days_list, self.__inf_list, label = "Infected")
        plt.plot(self.__days_list, self.__rec_list, label = "Recovered")
        ax.set_ylabel("Fraction of population")
        ax.set_xlabel("Days")
        plt.legend()
        plt.show()

