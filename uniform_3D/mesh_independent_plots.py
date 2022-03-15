from signac import get_project
from firedrake import *
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import re

plt.rcParams.update({
    "text.usetex": False,
    "font.family": "sans-serif",
    "font.sans-serif": ["Corbel"],
    "font.size" : 20
})

project = get_project("./")

#from init import nrefs

jobs_balance = list(project.find_jobs({'nref': {'$lt' : 5}}))
jobs_balance.sort(key=lambda job: job.sp.nref, reverse=True)
_, ax = plt.subplots()
#ax.set_prop_cycle(color=sns.color_palette("coolwarm_r",len(nrefs)))
colors = ['lime', 'green', 'springgreen', 'aquamarine']

for job, color in zip(jobs_balance, colors[0:len(jobs_balance)]):
    nref = job.sp.nref

    def read_output(output_file):
        obj_iteration_history_current = []
        iter_history_current = []

        with open(output_file, "r") as L2_file:
            for line in L2_file:
                obj_iteration = re.findall("obj: (-?[0-9]\d*)(\.\d+)?(e-?\+?\d*)?", line)
                current_it = re.findall("^It: (\d+)", line)
                if current_it:
                    iter_history_current.append(int(current_it[0]))
                if obj_iteration:
                    obj_it_float = float(obj_iteration[0][0] + obj_iteration[0][1] + obj_iteration[0][2])
                    obj_iteration_history_current.append(obj_it_float)

        return obj_iteration_history_current, iter_history_current

    L2_file = f"{job.ws}/output_{job.sp.nref}_L2.txt"
    obj_history, iter_history = read_output(L2_file)
    plt.loglog(iter_history, obj_history, '--', color=color, label=f"{nref}")

    euclidean_file = f"{job.ws}/output_{job.sp.nref}_euclidean.txt"
    obj_history, iter_history = read_output(euclidean_file)
    plt.loglog(iter_history, obj_history, '-', color=color, label=f"{nref}")


ax.set_xlabel('Number of iterations', fontsize=12)
ax.set_ylabel('Cost function', fontsize=12)

plt.setp(ax.get_xticklabels(), usetex=True)
plt.setp(ax.get_yticklabels(), usetex=True)
plt.grid(True)
plt.legend(loc="upper left", bbox_to_anchor=(1.05, 1))
plt.show()
# plt.savefig(
#     f"mesh_independent_plot.png",
#     dpi=1600,
#     orientation="portrait",
#     papertype=None,
#     format=None,
#     transparent=True,
#     bbox_inches="tight",
# )
