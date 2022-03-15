# project.py
import flow
from flow import FlowProject
import itertools
import platform
import sys
sys.settrace

plt = platform.system()


inner_products = ["L2", "euclidean"]

@FlowProject.label
def check_finished(job):
    check = True
    ref_str = job.sp.ref_strategy
    for inner_product in inner_products:
        check = check and job.isfile(f"finished_{ref_str}_{inner_product}_{job.sp['coarsen_refine']}.txt") \
                and job.isfile(f"output_{ref_str}_{inner_product}_{job.sp['coarsen_refine']}.txt")
    return check


@FlowProject.operation
@flow.cmd
@FlowProject.post(check_finished)
def launch_opti(job):
    if plt == "Linux":
        simulation = "source /g/g92/miguel/workspace/firedrake_latest.sh \n"
    else:
        simulation = ""

    ref_str = job.sp.ref_strategy

    for inner_product in inner_products:
        output =  f"{job.ws}/output_{ref_str}_{inner_product}_{job.sp['coarsen_refine']}.txt"
        if plt == "Linux":
            simulation += f"srun --output={output} "
        simulation += f"python3 {job.sp.case}.py \
                --output_dir {job.ws} \
                --ref_strategy {ref_str}\
                --inner_product {inner_product}\
                --coarsen_refine {job.sp['coarsen_refine']}"
        if plt == "Linux":
            simulation += "\n"
        else:
            simulation += f" &> {job.ws}/output_{ref_str}_{inner_product}_{job.sp['coarsen_refine']}.txt\n"
    return simulation


@FlowProject.label
def check_plots(job):
    return job.isfile("convergence_plots_1.svg") and job.isfile("convergence_plots_2.svg") \
           and job.isfile("convergence_plots_3.svg")


@FlowProject.operation
@flow.cmd
@FlowProject.pre(check_finished)
@FlowProject.post(check_plots)
def plot_convergence(job):
    post_process = "srun " if plt == "Linux" else ""
    cost_log = "--cost_log" if job.sp['case'] != 'thermal_flow' else ""
    post_process += f"python3 plot_history.py \
            {cost_log} \
            --dir {job.ws} \
            --ref_strategy {job.sp.ref_strategy}\
            --coarsen_refine {job.sp['coarsen_refine']}"
    return post_process


def filename_job(job, ip):
    return f"design_{job.sp.ref_strategy}_{ip}_{job.sp['coarsen_refine']}"


def screenshot_name(job, ip):
    return f"{job.sp.case}_{filename_job(job, ip)}_{job.sp.movlim}"


@FlowProject.label
def check_figures(job):
    check = True
    for inner_product in inner_products:
        check = check and job.isfile(f"{screenshot_name(job, inner_product)}.png") \
                    and job.isfile(f"{screenshot_name(job, inner_product)}_grid.png")
    return check


@FlowProject.operation
@flow.cmd
@FlowProject.pre(check_finished)
@FlowProject.post(check_figures)
def plot_figures(job):
    post_process = "srun " if plt == "Linux" else ""
    for inner_product in inner_products:
        post_process += f"pvpython screenshot_design.py \
                        --filename {filename_job(job, inner_product)} \
                        --output_name {screenshot_name(job, inner_product)} \
                        --results_dir {job.ws}\n"

    return post_process


if __name__ == '__main__':
    FlowProject().main()
