# project.py
import flow
from flow import FlowProject
import itertools
import platform
import sys
sys.settrace

plt = platform.system()
meshes = [0, 1]
inner_products = ["L2", "euclidean"]

@FlowProject.label
def check_finished(job):
    check = True
    for mesh, inner_product in itertools.product(
        meshes, inner_products
    ):
        check = check and job.isfile(f"finished_{mesh}_{inner_product}.txt") and job.isfile(f"output_{mesh}_{inner_product}.txt")
    return check


@FlowProject.operation
@flow.cmd
@FlowProject.post(check_finished)
def launch_opti(job):
    ref_dict = {"compliance": {0: 0, 1: 4},
                "stress": {0: 0, 1: 4},
                "mechanism": {0: 0, 1: 2}}
    simulation = ""

    for mesh, inner_product in itertools.product(
        meshes, inner_products
    ):
        output =  f"{job.ws}/output_{mesh}_{inner_product}.txt"
        simulation += f"python3 {job.sp.case}.py \
                --output_dir {job.ws} \
                --uniform {mesh}\
                --inner_product {inner_product} \
                --nref {ref_dict[job.sp.case][mesh]} "
        if job.sp.case == "stress":
            simulation += f" --quad_degree {job.sp.quad_degree} \
                            --p_norm_coeff {job.sp.p_norm_coeff} \
                            --penalty_param {job.sp.penalty_param} "
        simulation += f" > {job.ws}/output_{mesh}_{inner_product}.txt\n"
    return simulation


@FlowProject.label
def check_plots(job):
    return job.isfile("convergence_plots_1.svg") and job.isfile("convergence_plots_2.svg") and job.isfile("convergence_plots_3.svg")


@FlowProject.operation
@flow.cmd
@FlowProject.pre(check_finished)
@FlowProject.post(check_plots)
def plot_convergence(job):
    post_process = ""

    cost_log = "--cost_log" if job.sp['case'] != 'mechanism' else ""
    post_process += f"python3 plot_history.py \
            {cost_log} \
            --dir {job.ws}"
    return post_process


@FlowProject.label
def check_figures(job):
    check = True
    for mesh, inner_product in itertools.product(
        meshes, inner_products
    ):
        check = check and job.isfile(f"design_{mesh}_{inner_product}.png")
    return check


@FlowProject.operation
@flow.cmd
@FlowProject.pre(check_finished)
@FlowProject.post(check_figures)
def plot_figures(job):
    post_process = ""

    for mesh, inner_product in itertools.product(
        meshes, inner_products
    ):
        post_process += f"pvpython screenshot_design.py \
                --filename design_{mesh}_{inner_product} \
                --results_dir {job.ws}\n"

    return post_process


if __name__ == '__main__':
    FlowProject().main()
