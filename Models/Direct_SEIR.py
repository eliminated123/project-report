import matplotlib.pyplot as plt

class SEIRModel:
    
    def __init__(self):
        self.exposure_rate = 1
        self.incubation_rate = 1
        self.recovery_rate = 0.1
        self.__days_list = range(0,101)
        self.__sus_list = [0.99]
        self.__exp_list = [0.01]
        self.__inf_list = [0]
        self.__rec_list = [0]

    def __update_sus (self, sus, inf):
        return sus - self.exposure_rate * inf * sus

    def __update_exp (self, exp, sus, inf):
        return exp + self.exposure_rate * inf * sus - self.incubation_rate * exp 

    def __update_inf (self, inf, exp):
        return inf + self.incubation_rate * exp - self.recovery_rate * inf

    def __update_rec (self, rec, inf):
        return rec + self.recovery_rate * inf

    def calculate(self):
        for i in range(0, 100):
            sus = self.__sus_list[-1]
            exp = self.__exp_list[-1]
            inf = self.__inf_list[-1]
            rec = self.__rec_list[-1]
            
            self.__sus_list.append(self.__update_sus(sus, inf))
            self.__exp_list.append(self.__update_exp(exp, sus, inf))
            self.__inf_list.append(self.__update_inf(inf, exp))
            self.__rec_list.append(self.__update_rec(rec, inf))

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

model = SEIRModel()
model.calculate()
model.plot()