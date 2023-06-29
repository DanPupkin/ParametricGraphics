import numpy as np
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import sympy as smp
from sympy.abc import t
from sympy import Interval, oo, im, periodicity, is_monotonic, singularities, is_decreasing, is_increasing

class basic_visualisator:
    def __init__(self, start=-5, end=5, accuracy: int = 300):
        self.start = start
        self.end = end
        self.t_range = np.linspace(start, end, accuracy, dtype=np.float32)
        self.accuracy = accuracy
        self.epsilon = 0.001
        self.x_view_start, self.x_view_end = -3, 3
        self.y_view_start, self.y_view_end = -3, 3
    def two_var_parametric_function_print(self, first_function, second_function, argument = "t"):
        first_excr = smp.sympify(first_function)
        second_excr = smp.sympify(second_function)
        period = periodicity(first_excr, t, check=True)
        if period:
            print(f"функция переодическая с периодом {period}")
            if abs(self.end - self.start) > period:
                self.end = self.start + 2*period
#из-за симметричности не всегда хватает прибавить 1 период, требуется 2, а вводить проверку на симметричность не представляется
#возможным

        xt_breaks = singularities(first_excr, t)
        yt_breaks = singularities(second_excr, t)
        for i in range(len(self.t_range)):
            if self.t_range[i] in xt_breaks:
                self.t_range[i] += self.epsilon
        xt = np.array([first_excr.subs(argument, a) for a in self.t_range], dtype=np.float32)
        yt = np.array([second_excr.subs(argument, a) for a in self.t_range], dtype=np.float32)
        x_view_start, x_view_end = self.x_view_start, self.x_view_end
        y_view_start, y_view_end = self.y_view_start, self.y_view_end
        finally_plots = []
        xt_solves = []
        #производная x по t и ее нули
        xt_diff = smp.diff(first_excr, t)
        yt_diff = smp.diff(second_excr, t)
        # вычисляем 1 и 2 производну по x
        yx_first_diff = yt_diff / xt_diff
        yx_first_diff = yx_first_diff.simplify()
        yx_second_diff = smp.diff(yx_first_diff, t) / xt_diff
        yx_second_diff = yx_second_diff.simplify()
        yx_second_diff_breaks = singularities(yx_second_diff, t, domain=Interval(self.start, self.end))
        yx_second_diff_solves = smp.solveset(yx_second_diff, domain=Interval(self.start, self.end))
        inflection_points = yx_second_diff_breaks | yx_second_diff_solves
        print(yx_second_diff.simplify(), "yx_second_diff")
        #xt_diff_solves = smp.solve(xt_diff, t)
        xt_diff_solves = smp.solveset(xt_diff, domain = Interval(self.start, self.end))
        #сохраняем все особые точки и нули производной
        xt_diff_singularities = singularities(xt_diff, t)
        print(xt_diff_singularities, "xt_diff singularities")
        xt_monotonic_breakpoints = set(xt_diff_singularities)
        [xt_monotonic_breakpoints.add(a) for a in xt_breaks]
        [xt_monotonic_breakpoints.add(a) for a in xt_diff_solves]
        print(xt_diff_solves, "diff solves")
        xt_monotonic_breakpoints.add(self.start)
        xt_monotonic_breakpoints.add(self.end)
        #чистка от комплексных чисел (ноль убрали, потом вернули)
        if 0 in xt_monotonic_breakpoints:
            xt_monotonic_breakpoints.remove(0)
            xt_monotonic_breakpoints = set([i for i in xt_monotonic_breakpoints if im(i)==0] )
            xt_monotonic_breakpoints.add(0)
        else:
            xt_monotonic_breakpoints = set([i for i in xt_monotonic_breakpoints if im(i)==0])
        xt_monotonic_breakpoints = sorted(xt_monotonic_breakpoints)

        xt_monotonic_intervals = []
        #создаем интервалы монотонности
        for point_int in range(len(xt_monotonic_breakpoints)-1):
            if xt_monotonic_breakpoints[point_int+1]>=self.start and xt_monotonic_breakpoints[point_int+1]<=self.end:
                # (a; b)
                if is_monotonic(first_excr, Interval(xt_monotonic_breakpoints[point_int], xt_monotonic_breakpoints[point_int+1])):
                    xt_monotonic_intervals.append(Interval(xt_monotonic_breakpoints[point_int], xt_monotonic_breakpoints[point_int+1]))
                    #Lopen= [a; b)
                elif is_monotonic(first_excr, Interval.Lopen(xt_monotonic_breakpoints[point_int], xt_monotonic_breakpoints[point_int+1])):
                    xt_monotonic_intervals.append(Interval.Lopen(xt_monotonic_breakpoints[point_int], xt_monotonic_breakpoints[point_int+1]))
                    #Ropen= (a; b]
                elif is_monotonic(first_excr, Interval.Ropen(xt_monotonic_breakpoints[point_int], xt_monotonic_breakpoints[point_int+1])):
                    xt_monotonic_intervals.append(Interval.Ropen(xt_monotonic_breakpoints[point_int], xt_monotonic_breakpoints[point_int+1]))
                    #open= [a;b]
                elif is_monotonic(first_excr, Interval.open(xt_monotonic_breakpoints[point_int], xt_monotonic_breakpoints[point_int+1])):
                    xt_monotonic_intervals.append(Interval.open(xt_monotonic_breakpoints[point_int], xt_monotonic_breakpoints[point_int+1]))

        #interval.args[0] and[1] is a start and end of the interval
        print(xt_monotonic_intervals, "intervals")

        list_of_asymptotes = []
        for interval in xt_monotonic_intervals:
            print(interval)
            left_border = interval.args[0]
            right_border = interval.args[1]
            left_bracket = "[" if interval.args[2] else "("
            right_bracket = "]" if interval.args[3] else ")"
            sub_xt = np.array([first_excr.subs(argument, a) for a in self.t_range if a < right_border and a > left_border])
            sub_yt = np.array(
                [second_excr.subs(argument, a) for a in self.t_range if a < right_border and a > left_border])
            for i in range(len(sub_xt)):
                    for brk in yt_breaks:
                        if sub_xt[i] != None:
                            if abs(sub_xt[i]-brk) <=0.04:
                                sub_xt[i] = None


            fig = plt.figure()
            ax = fig.add_subplot()
            fig.suptitle(f'Ветвь кривой для t ∈ {left_bracket}{left_border};{right_border}{right_bracket}', fontsize=14, fontweight='bold')
            finally_plots.append((sub_xt, sub_yt))
            plt.plot(sub_xt, sub_yt)
            x_label = ""
            if interval.args[0] in xt_breaks:
                k_koeff = smp.limit(second_excr/first_excr, t, left_border, "+")
                if k_koeff != oo and k_koeff != -oo:
                    b_koeff = smp.limit(second_excr - k_koeff*first_excr, t, left_border, "+")
                    asymptote = t*k_koeff + b_koeff
                    asymptote_x = np.linspace(-100, 100, 5)
                    asymptote_y = np.array([asymptote.subs(t, a) for a in asymptote_x])
                    list_of_asymptotes.append((asymptote_x, asymptote_y))
                    x_label += f'имеется асимптота {k_koeff}x + ({b_koeff}) при t→{left_border}+0\n'
                    ax.set_xlabel(x_label)
                    plt.plot(asymptote_x, asymptote_y, linestyle="--")

            if interval.args[1] in xt_breaks:
                k_koeff = smp.limit(second_excr/first_excr, t, right_border, "-")
                if k_koeff != oo and k_koeff != -oo:
                    b_koeff = smp.limit(second_excr - k_koeff*first_excr, t, right_border, "-")
                    asymptote = t*k_koeff + b_koeff
                    asymptote_x = np.linspace(-100, 100, 5)
                    asymptote_y = np.array([asymptote.subs(t, a) for a in asymptote_x])
                    list_of_asymptotes.append((asymptote_x, asymptote_y))
                    x_label += f'имеется асимптота {k_koeff}x + ({b_koeff}) при t→{right_border}-0\n'
                    ax.set_xlabel(x_label)
                    plt.plot(asymptote_x, asymptote_y, linestyle="--")

            plt.ylim(y_view_start, y_view_end)
            plt.xlim(x_view_start, x_view_end)

            #print(smp.solveset(yx_second_diff, domain=Interval.open(self.start, self.end)), "sec diff solves")
            a = is_decreasing(yx_first_diff, interval)
            b = is_increasing(yx_first_diff, interval)
            if is_monotonic(yx_first_diff, interval) or (bool(a) ^ bool(b)):
                #монотонность позволяет проверить знак проверив середину интервала
                znak = ">" if yx_second_diff.subs(t, left_border+(right_border-left_border)/2) >0 else "<"
                if znak == ">":
                    x_label += f'вторая производная {znak} 0, функция выпукла вниз\n'
                else:
                    x_label += f'вторая производная {znak} 0, функция выпукла вверх\n'
                ax.set_xlabel(x_label)
            ax = plt.gca()
            ax.axhline(y=0, color='k')
            ax.axvline(x=0, color='k')
            plt.show()
        fig = plt.figure()
        for plott in finally_plots:
            plot_x, plot_y = plott
            plt.plot(plot_x, plot_y, color="blue")
        for asymp in list_of_asymptotes:
            plt.plot(asymp[0],asymp[1], linestyle='--')
        for point in inflection_points:
            if point not in xt_breaks and point not in yt_breaks:
                        #перегиб можно пометить стрелкой или точкой, а можно и так и так
                        #чтобы не перегружать рисунки пометим перегибы цветом
                        #plt.annotate("Точка перегиба", xy=(first_excr.subs(t, point), second_excr.subs(t, point)), xytext=(first_excr.subs(t, point), second_excr.subs(t, point)+2), arrowprops=dict(arrowstyle='->'))
                plt.scatter(first_excr.subs(t, point), second_excr.subs(t, point), color="red", marker="*", s=100)
                print(first_excr.subs(t, point), second_excr.subs(t, point), "точка перегиба")
        plt.ylim(y_view_start, y_view_end)
        plt.xlim(x_view_start, x_view_end)
        ax = plt.gca()
        ax.axhline(y=0, color='k')
        ax.axvline(x=0, color='k')
        blue_path = mpatches.Patch(color='blue', label='график y(x)')
        dotted_line =  mlines.Line2D([], [], color='black', linestyle="--", label='асимптоты')
        inflect_mark = mlines.Line2D([0], [0], color='red', marker="*", label='точки перегиба')
        plt.legend(handles=[blue_path, dotted_line, inflect_mark])
        fig.suptitle(f"x(t)={first_function}, y(t)={second_function}", fontsize=14, fontweight='bold')
        plt.show()


def main():
    #невозможно считать значения на (-inf;inf)
    #из-за чего на интервалах, где функция xt монотонна могут быть разрывы, но с увелечением диапазона t
    #разрывы становятся меньше или пропадают
    visual = basic_visualisator(-16, 16, 300)

    #x = "(t**2)/(t**2-1)" Пример 19
    #y = "(t**2)/(t**3+1)"
    x = "t/(1+t**3)" #Декартов лист
    y = "(t**2)/(1+t**3)"
    list_of_functions=[
        ("t/(t**2-1)", "(t**2)/(t-1)"),                         #0-367
        ("(4-t**2)/(1+t**3)","(t**2)/(1+t**3)"),                #1-368
        ("(t**3)/(t**3+1)","(t**2)/(1+t**3)"),                  #2-369
        ("t**2","(t**3+2t**2+t)/(t+2)"),                        #3-370
        ("(t**2-1)/(t*(t+2))","(t**2)/((t+2)*(t+1))"),          #4-371
        ("(t**2+1)/(4*(t-1))","t/(t+1)"),                       #5-372
        ("t/(1-t**2)","(t*(1-4*t**2))/(1-t**2)"),               #6-373
        ("(t**2)/(t**2-1)","(t**2+1)/(t+2)")                    #7-374
    ]
    #x,y = list_of_functions[0]
    visual.two_var_parametric_function_print(x, y)

if __name__ == '__main__':
    main()
