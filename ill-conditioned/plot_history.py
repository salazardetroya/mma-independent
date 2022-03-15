import re
import sys
import matplotlib.pyplot as plt

plt.rcParams["text.usetex"] = True
import numpy as np
import glob
import re
import argparse

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument(
    "--cost_log",
    action="store_true",
    dest="cost_log",
    help="Plot cost function in loglog",
    default=False,
)
parser.add_argument(
    "--dir",
    action="store",
    dest="folder_dir",
    type=str,
    help="Folder of results to plot",
    default="./",
)
parser.add_argument(
    "--tikz",
    action="store",
    dest="tikz",
    type=bool,
    help="Print to tikz figures",
    default=False,
)
parser.add_argument(
    "--iters",
    action="store",
    dest="iters",
    type=int,
    help="Max iters to plot",
    default=200,
)
args = parser.parse_args()
cost_log = args.cost_log
tikz = args.tikz
folder_dir = args.folder_dir
iters = args.iters

files = glob.glob(folder_dir + "/*.txt")

needed_files = [
    "output_1_L2.txt",
    "output_1_euclidean.txt",
    "output_0_L2.txt",
    "output_0_euclidean.txt",
]
labels = [
    "Uniform in $L^2$",
    "Uniform in $\ell^2$",
    "Non-uniform in $L^2$",
    "Non-uniform in $\ell^2$",
]

result_files = []
for needed_file in needed_files:
    regex = re.compile(".*" + needed_file)
    file_found = False
    for s in files:
        reg_result = regex.search(s)
        if reg_result:
            result_files.append(reg_result.group(0))
            file_found = True
            break
    assert file_found, "We need to have file " + needed_file

uniform_true_file = result_files[0]
uniform_false_file = result_files[1]
high_true_file = result_files[2]
high_false_file = result_files[3]

# There should be a way to get this all together in a single loop
constr_history = []
constr2_history = []
obj_history = []
iter_history = []
kkt_conditions = []
for file_output in [
    uniform_true_file,
    uniform_false_file,
    high_true_file,
    high_false_file,
]:
    obj_iteration_history = []
    iter_history_current = []
    constr_iteration_history = []
    constr2_iteration_history = []
    kkt_history = []

    file_object = open(file_output, "r")
    current_it = 0
    for line in file_object:
        reg_number = "(-?[0-9]\d*)(\.\d+)?(e-?\+?\d*)?"
        obj_iteration = re.findall(f"obj: {reg_number}", line)
        constraint_itera = re.findall("g\[0\]: (-?[0-9]\d*)(\.\d+)?(e-?\+?\d*)?", line)
        constraint_itera_2 = re.findall(
            "g\[1\]: (-?[0-9]\d*)(\.\d+)?(e-?\+?\d*)?", line
        )
        kkt_iter = re.findall("kkt: (-?[0-9]\d*)(\.\d+)?(e-?\+?\d*)?", line)
        if obj_iteration:
            current_it += 1
            iter_history_current.append(current_it)
            obj_it_float = float(
                obj_iteration[0][0] + obj_iteration[0][1] + obj_iteration[0][2]
            )
            obj_iteration_history.append(obj_it_float)
        if constraint_itera:
            constr_itera = float(
                constraint_itera[0][0] + constraint_itera[0][1] + constraint_itera[0][2]
            )
            constr_iteration_history.append(constr_itera)
        if constraint_itera_2:
            constr_itera = float(
                constraint_itera_2[0][0]
                + constraint_itera_2[0][1]
                + constraint_itera_2[0][2]
            )
            constr2_iteration_history.append(constr_itera)
        if kkt_iter:
            kkt_it_float = float(kkt_iter[0][0] + kkt_iter[0][1] + kkt_iter[0][2])
            kkt_history.append(float(kkt_it_float))

    obj_history.append(np.array(obj_iteration_history))
    constr_history.append(np.array(constr_iteration_history))
    kkt_conditions.append(np.array(kkt_history))
    if constr2_iteration_history:
        constr2_history.append(np.array(constr2_iteration_history))
    iter_history.append(iter_history_current)

## Remove repeated indices from having restarted a simulation
# for i in range(4):
#    iter_history[i], indices = np.unique(iter_history[i], return_index=True)
#    obj_history[i] = obj_history[i][indices]
#    constr_history[i] = constr_history[i][indices]
#    if constr2_history:
#        constr2_history[i] = constr2_history[i][indices]
#    # kkt_conditions[i] = kkt_conditions[i][indices]


fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
fig4, ax4 = plt.subplots()
if constr2_history:
    fig3, ax3 = plt.subplots()
figures_list = [fig1, fig2, fig4]
if constr2_history:
    figures_list.append(fig3)


color_to_plot = ["b-", "b--", "r-", "r--"]
labels = [
    "Uniform in $L^2$",
    "Uniform in $\ell^2$",
    "Non-uniform in $L^2$",
    "Non-uniform in $\ell^2$",
]

iter_cut = iters
for i, my_color, label in zip(range(4), color_to_plot, labels):
    # kkt_convergence_index = np.argmax(kkt_conditions[i][iter_cut:] < 1e-3)
    kkt_convergence_index = 1
    print("Last iteration to plot is: {}".format(iter_cut + kkt_convergence_index))
    if kkt_convergence_index > 0:
        idxs = list(range(iter_cut + kkt_convergence_index))
    else:
        idxs = list(range(800))

    if cost_log:
        ax1.loglog(
            iter_history[i][idxs[0] : idxs[-1]],
            obj_history[i][idxs[0] : idxs[-1]],
            my_color,
            label=label,
        )
    else:
        ax1.plot(
            iter_history[i][idxs[0] : idxs[-1]],
            obj_history[i][idxs[0] : idxs[-1]],
            my_color,
            label=label,
        )

    ax2.plot(
        iter_history[i][idxs[0] : idxs[-1]],
        constr_history[i][idxs[0] : idxs[-1]],
        my_color,
        label=label,
    )
    if constr2_history:
        ax3.plot(
            iter_history[i][idxs[0] : idxs[-1]],
            constr2_history[i][idxs[0] : idxs[-1]],
            my_color,
            label=label,
        )
    ax4.semilogy(
        iter_history[i][idxs[0] : idxs[-1]],
        kkt_conditions[i][idxs[0] : idxs[-1]],
        my_color,
        label=label,
    )


ax1.set_xlabel("Number of iterations", fontsize=18)
ax1.set_ylabel("Cost function", fontsize=18)
legend = ax1.legend(loc="best")
for label in legend.get_texts():
    label.set_fontsize(18)

ax2.set_xlabel("Number of iterations", fontsize=18)
ax2.set_ylabel("Constraint", fontsize=18)
legend = ax2.legend(loc="best")
for label in legend.get_texts():
    label.set_fontsize(18)

ax4.set_xlabel("Number of iterations", fontsize=18)
ax4.set_ylabel("Stopping criteria", fontsize=18)
legend = ax4.legend(loc="best")
for label in legend.get_texts():
    label.set_fontsize(18)

if constr2_history:
    ax3.set_xlabel("Number of iterations", fontsize=18)
    ax3.set_ylabel("Constraint", fontsize=18)
    legend = ax3.legend(loc="best")
    for label in legend.get_texts():
        label.set_fontsize(18)


if tikz:
    from tikzplotlib import save as tikz_save

    for figura in figures_list:
        plt.figure(figura.number)
        tikz_save(
            folder_dir + "/convergence_plots{}".format(figura.number) + ".tex",
            axis_height="\\figureheight",
            axis_width="\\figurewidth",
        )
else:
    for figura in figures_list:
        plt.figure(figura.number)
        plt.savefig(
            f"{folder_dir}/convergence_plots_{figura.number}.svg",
            dpi=1600,
            orientation="portrait",
            papertype=None,
            format=None,
            transparent=True,
            bbox_inches="tight",
        )

