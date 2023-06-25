import numpy as np
import matplotlib.pyplot as plt
import sympy as smp
from sympy.abc import t, x
from sympy import is_monotonic
from sympy import singularities
from sympy import Interval, oo
from sympy import im
class basic_visualisator:
    def __init__(self, start, end, accuracy:int):
        self.start = start
        self.end = end
        self.t_range = np.linspace(start, end, accuracy, dtype=np.float32)
        self.accuracy = accuracy
        self.epsilon = 0.001
    def function_print(self, func:str, argument:str, point=2):
        func_excr = smp.sympify(func)
        print(f"func:{func_excr}, point: {point}, value={func_excr.subs(argument, point)}")
    def two_var_parametric_function_print(self, first_function, second_function, argument = "t"):
        first_excr = smp.sympify(first_function)
        second_excr = smp.sympify(second_function)
        xt_breaks = singularities(first_excr, t)
        for i in range(len(self.t_range)):
            if self.t_range[i] in xt_breaks:
                self.t_range[i] += self.epsilon
        xt = np.array([first_excr.subs(argument, a) for a in self.t_range], dtype=np.float32)
        yt = np.array([second_excr.subs(argument, a) for a in self.t_range], dtype=np.float32)
        x_view_start, x_view_end = -5,5
        y_view_start, y_view_end = -5,5

        #точки разрыва и нули графика x от t

        #xt_solves = smp.solve(first_excr, t)
        xt_solves = []

        #производная x по t и ее нули
        xt_diff = smp.diff(first_excr, t)
        #xt_diff_solves = smp.solve(xt_diff, t)
        xt_diff_solves = smp.solveset(xt_diff, domain = Interval(self.start, self.end))
        #сохраняем все "особые" точки функции x от t в один лист
        xt_monotonic_breakpoints = set(xt_solves)
        [xt_monotonic_breakpoints.add(a) for a in xt_breaks]
        [xt_monotonic_breakpoints.add(a) for a in xt_diff_solves]
        print(xt_diff_solves, "diff solves")
        xt_monotonic_breakpoints.add(self.start)
        xt_monotonic_breakpoints.add(self.end)

        if 0 in xt_monotonic_breakpoints:
            xt_monotonic_breakpoints.remove(0)

            xt_monotonic_breakpoints = set([i for i in xt_monotonic_breakpoints if im(i)==0] )
            xt_monotonic_breakpoints.add(0)
        else:
            xt_monotonic_breakpoints = set([i for i in xt_monotonic_breakpoints if im(i)==0])
        xt_monotonic_breakpoints = sorted(xt_monotonic_breakpoints)

        xt_monotonic_intervals = []

        for point_int in range(len(xt_monotonic_breakpoints)-1):
            if xt_monotonic_breakpoints[point_int+1]>=self.start and xt_monotonic_breakpoints[point_int+1]<=self.end:
                # (a; b)
                if is_monotonic(first_excr, Interval(xt_monotonic_breakpoints[point_int], xt_monotonic_breakpoints[point_int+1])):
                    xt_monotonic_intervals.append(Interval(xt_monotonic_breakpoints[point_int], xt_monotonic_breakpoints[point_int+1]))
                    #Lopen [a; b)
                elif is_monotonic(first_excr, Interval.Lopen(xt_monotonic_breakpoints[point_int], xt_monotonic_breakpoints[point_int+1])):
                    xt_monotonic_intervals.append(Interval.Lopen(xt_monotonic_breakpoints[point_int], xt_monotonic_breakpoints[point_int+1]))
                    #Ropen (a; b]
                elif is_monotonic(first_excr, Interval.Ropen(xt_monotonic_breakpoints[point_int], xt_monotonic_breakpoints[point_int+1])):
                    xt_monotonic_intervals.append(Interval.Ropen(xt_monotonic_breakpoints[point_int], xt_monotonic_breakpoints[point_int+1]))
                    #open [a;b]
                elif is_monotonic(first_excr, Interval.open(xt_monotonic_breakpoints[point_int], xt_monotonic_breakpoints[point_int+1])):
                    xt_monotonic_intervals.append(Interval.open(xt_monotonic_breakpoints[point_int], xt_monotonic_breakpoints[point_int+1]))

        #interval.args[0] and[1] is a start and end of the interval
        print(xt_monotonic_intervals, "intervals")

# нормализация графика на разрывах
        #так же где-то в этом районе надо будет позднее дописать фикс того, что иксы пропадают, когда
        #их очень много в окресностях одной точки
        #xt[xt.argmax()] = None
        #xt[xt.argmin()] = None

        for i in range(len(xt)-1):
            if xt[i]:
                if abs(xt[i] - xt[i+1])>= 5:
                    xt[i] = None


        for interval in xt_monotonic_intervals:
            print(interval)
            left_border = interval.args[0]
            right_border = interval.args[1]
            sub_xt = np.array([first_excr.subs(argument, a) for a in self.t_range if a < right_border and a > left_border])
            sub_yt = np.array([second_excr.subs(argument, a) for a in self.t_range if a < right_border and a > left_border])
            plt.plot(sub_xt, sub_yt)

            if interval.args[0] in xt_breaks:
                k_koeff = smp.limit(second_excr/first_excr, t, interval.args[0], "+")
                if k_koeff != oo and k_koeff != -oo:
                    b_koeff = smp.limit(second_excr - k_koeff*first_excr, t, interval.args[0], "+")
                    kasatelnaya = t*k_koeff + b_koeff
                    print(interval.args)
                    kasatelnaya_x = np.linspace(-100, 100, 10)
                    kasatelnaya_y = np.array([kasatelnaya.subs(t, a) for a in kasatelnaya_x])
                    plt.plot(kasatelnaya_x, kasatelnaya_y)

            if interval.args[1] in xt_breaks:
                k_koeff = smp.limit(second_excr/first_excr, t, interval.args[1], "-")
                if k_koeff != oo and k_koeff != -oo:
                    b_koeff = smp.limit(second_excr - k_koeff*first_excr, t, interval.args[1], "-")
                    kasatelnaya = t*k_koeff + b_koeff
                    print(interval.args)
                    kasatelnaya_x = np.linspace(-100, 100, 10)
                    kasatelnaya_y = np.array([kasatelnaya.subs(t, a) for a in kasatelnaya_x])
                    plt.plot(kasatelnaya_x, kasatelnaya_y)

            plt.ylim(y_view_start, y_view_end )
            plt.xlim(x_view_start, x_view_end )

            ax = plt.gca()
            ax.axhline(y=0, color='k')
            ax.axvline(x=0, color='k')
            plt.show()



        plt.plot(xt, yt)
        plt.ylim(y_view_start, y_view_end)
        plt.xlim(x_view_start, x_view_end)
        ax = plt.gca()
        ax.axhline(y=0, color='k')
        ax.axvline(x=0, color='k')
        plt.show()


def main():
    #невозможно считать значения на (-inf;inf), поэтому нужно как-то определить оптимальную длину интервала, которой хватит всем
    visual = basic_visualisator(-10, 10, 300)

    #x = "(t**2)/(t**2-1)"
    #y = "(t**2)/(t**3+1)"
    #x = "t/(1+t**3)"
    #y = "(t**2)/(1+t**3)"
    x = "2*sin(2*t)"
    y = "2*cos(t)"
    visual.two_var_parametric_function_print(x, y)

if __name__ == '__main__':
    main()
