import KratosMultiphysics
import matplotlib.pyplot as plt
import numpy as np

class PlotManufacturedFieldUtility(object):
    def __init__(self, settings = KratosMultiphysics.Parameters("""{}""")):

        default_settings = KratosMultiphysics.Parameters("""{
            "benchmark_name"  : "custom_body_force.polynomial_vortex",
            "x_range"         : [0.0, 1.0],
            "y_range"         : [0.0, 1.0],
            "Parameters"      : {}
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        # from polynomial_vortex import ManufacturedSolution
        benchmark_module = __import__(settings["benchmark_name"].GetString(), fromlist=[None])
        self.benchmark = benchmark_module.CreateManufacturedSolution(settings["Parameters"])

        L = 0.5
        n_points = 30
        step = L / n_points
        x = np.arange(settings["x_range"][0].GetDouble(), settings["x_range"][1].GetDouble(), step)
        y = np.arange(settings["y_range"][0].GetDouble(), settings["y_range"][1].GetDouble(), step)
        self.X, self.Y = np.meshgrid(x, y)

    def PlotVelocity(self, time = 1.0):
        fig, ax = plt.subplots()
        self.plot_velocity(ax, time)
        plt.show()

    def PlotBodyForce(self, time = 1.0):
        fig, ax = plt.subplots()
        self.plot_body_force(ax, time)
        plt.show()

    def DynamicPlotVelocity(self, start = 0.0, stop = 5.0, step = 0.1, pause = 0.5):
        fig, ax = plt.subplots()
        for time in np.arange(start, stop, step):
            self.plot_velocity(ax, time)
            plt.pause(pause)

    def plot_velocity(self, ax, time):
        U1 = self.benchmark.u1(self.X, self.Y, time)
        U2 = self.benchmark.u2(self.X, self.Y, time)
        self.plot_vectors(ax, U1, U2)
    
    def plot_body_force(self, ax, time):
        B1 = self.benchmark.body_force1(self.X, self.Y, time)
        B2 = self.benchmark.body_force2(self.X, self.Y, time)
        self.plot_vectors(ax, B1, B2)

    def plot_vectors(self, ax, U, V):
        ax.cla()
        q = ax.quiver(self.X, self.Y, U, V)
        ax.quiverkey(q, X=0.3, Y=1.1, U=1, label='Quiver key, length = 1', labelpos='E')

if __name__ == "__main__":
    PlotManufacturedFieldUtility().PlotVelocity()
